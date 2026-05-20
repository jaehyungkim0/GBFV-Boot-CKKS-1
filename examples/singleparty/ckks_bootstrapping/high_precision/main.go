// GBFV Bootstrapping from CKKS, relying on META-BTS to instantiate high precision CKKS bootstrapping.
//
// Use -short to run with smaller insecure parameters, and -once for a quick smoke test.
package main

import (
	"crypto/rand"
	"flag"
	"fmt"
	"math"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")
var flagOnce = flag.Bool("once", false, "run only one experiment iteration.")

const (
	defaultLogN           = 16
	shortLogNDelta        = 3
	defaultIterations     = 20
	gbfvBase              = int64(16)
	gbfvSlots             = int64(128)
	perturbationMagnitude = 1e-9
)

type bootstrappingSetup struct {
	params    ckks.Parameters
	btpParams bootstrapping.Parameters
}

type schemeContext struct {
	encoder   *ckks.Encoder
	encryptor *rlwe.Encryptor
	decryptor *rlwe.Decryptor
	evaluator *bootstrapping.Evaluator
}

type gbfvParameters struct {
	base            int64
	slots           int64
	ringDegree      int64
	baseInt         *big.Int
	upperBound      *big.Int
	upperBoundFloat *big.Float
	q               *big.Int
	qPrime          *big.Int
	qFloat          *big.Float
	qPrimeFloat     *big.Float
	oneFloat        *big.Float
	precision       uint
}

func main() {
	flag.Parse()

	logN := defaultLogN
	if *flagShort {
		logN -= shortLogNDelta
	}

	setup, err := newBootstrappingSetup(logN)
	if err != nil {
		panic(err)
	}

	printBootstrappingParameters(setup.btpParams)

	ctx, err := newSchemeContext(setup.params, setup.btpParams)
	if err != nil {
		panic(err)
	}

	gbfv := newGBFVParameters(setup.params, logN)
	printGBFVParameters(gbfv)

	runExperiment(setup.params, ctx, gbfv)
}

func newBootstrappingSetup(logN int) (bootstrappingSetup, error) {
	params, err := newResidualParameters(logN)
	if err != nil {
		return bootstrappingSetup{}, err
	}

	btpParams, err := newBootstrappingParameters(params, logN)
	if err != nil {
		return bootstrappingSetup{}, err
	}

	if *flagShort {
		// Corrects Q0/|m(X)| for the smaller number of slots while keeping the same precision target.
		btpParams.Mod1ParametersLiteral.LogMessageRatio += defaultLogN - params.LogN()
	}

	return bootstrappingSetup{params: params, btpParams: btpParams}, nil
}

func newResidualParameters(logN int) (ckks.Parameters, error) {
	return ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            []int{60, 45, 45, 45},
		LogP:            []int{61, 61, 61, 61},
		LogDefaultScale: 90,
		Xs:              ring.Ternary{H: 192},
	})
}

func newBootstrappingParameters(params ckks.Parameters, logN int) (bootstrapping.Parameters, error) {
	return bootstrapping.NewParametersFromLiteral(params, bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(logN),
		LogP: []int{61, 61, 61, 61},
		IterationsParameters: &bootstrapping.IterationsParameters{
			BootstrappingPrecision: []float64{25, 25, 5},
			ReservedPrimeBitSize:   28,
		},
		Xs: params.Xs(),
	})
}

func printBootstrappingParameters(btpParams bootstrapping.Parameters) {
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(),
		btpParams.ResidualParameters.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())
}

func newSchemeContext(params ckks.Parameters, btpParams bootstrapping.Parameters) (schemeContext, error) {
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		return schemeContext{}, err
	}
	fmt.Println("Done")

	evaluator, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		return schemeContext{}, err
	}

	return schemeContext{
		encoder:   encoder,
		encryptor: encryptor,
		decryptor: decryptor,
		evaluator: evaluator,
	}, nil
}

func newGBFVParameters(params ckks.Parameters, logN int) gbfvParameters {
	precision := params.EncodingPrecision()
	ringDegree := int64(1) << logN
	baseInt := big.NewInt(gbfvBase)

	upperBound := new(big.Int).Exp(baseInt, big.NewInt(ringDegree/gbfvSlots), nil)
	upperBound.Add(upperBound, big.NewInt(1))

	ringQ := params.RingQ().AtLevel(params.LevelsConsumedPerRescaling() - 1)
	q := new(big.Int).Set(ringQ.ModulusAtLevel[params.LevelsConsumedPerRescaling()-1])
	qPrime := new(big.Int).Div(new(big.Int).Set(ringQ.ModulusAtLevel[2*params.LevelsConsumedPerRescaling()-1]), q)

	return gbfvParameters{
		base:            gbfvBase,
		slots:           gbfvSlots,
		ringDegree:      ringDegree,
		baseInt:         baseInt,
		upperBound:      upperBound,
		upperBoundFloat: new(big.Float).SetPrec(precision).SetInt(upperBound),
		q:               q,
		qPrime:          qPrime,
		qFloat:          new(big.Float).SetPrec(precision).SetInt(q),
		qPrimeFloat:     new(big.Float).SetPrec(precision).SetInt(qPrime),
		oneFloat:        new(big.Float).SetPrec(precision).SetInt64(1),
		precision:       precision,
	}
}

func printGBFVParameters(gbfv gbfvParameters) {
	fmt.Println("Q = ", gbfv.q)
	fmt.Printf("GBFV Bootstrapping Parameters: b = %d, k = %d, log2(b^(N/k)) = %.2f\n",
		gbfv.base,
		gbfv.slots,
		plaintextPrecisionBits(gbfv))
}

func plaintextPrecisionBits(gbfv gbfvParameters) float64 {
	return float64(gbfv.ringDegree/gbfv.slots) * math.Log2(float64(gbfv.base))
}

func runExperiment(params ckks.Parameters, ctx schemeContext, gbfv gbfvParameters) {
	iterations := iterationCount()
	var totalDuration time.Duration

	for iter := 0; iter < iterations; iter++ {
		fmt.Println("-----------------------------------------------------")
		fmt.Printf("------------------ Iteration %d ----------------------\n", iter)
		fmt.Println("-----------------------------------------------------")

		duration := runIteration(params, ctx, gbfv, iter)
		totalDuration += duration
	}

	fmt.Println("-------------------------------------------------------------------")
	fmt.Println("Average Bootstrapping Time:", totalDuration/time.Duration(iterations))
}

func iterationCount() int {
	if *flagOnce {
		return 1
	}
	return defaultIterations
}

func runIteration(params ckks.Parameters, ctx schemeContext, gbfv gbfvParameters, iter int) time.Duration {
	randomInts, err := sampleRandomInts(gbfv)
	if err != nil {
		panic(err)
	}

	valuesScaled := gbfvCoefficientValues(randomInts, gbfv)
	valuesWant := normalizeByQ(valuesScaled, gbfv)
	valuesPerturbed := perturbValues(valuesScaled, gbfv)

	plaintext, ptT, ptInverseT, err := encodePlaintexts(params, ctx.encoder, gbfv, valuesPerturbed)
	if err != nil {
		panic(err)
	}

	ciphertext, err := ctx.encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	start := time.Now()
	ciphertextBefore, ciphertextAfter, err := bootstrapGBFV(ctx.evaluator, ciphertext, ptT, ptInverseT, gbfv)
	if err != nil {
		panic(err)
	}
	elapsed := time.Since(start)

	fmt.Printf("GBFV Bootstrapping Time (iteration %d): %s\n", iter, elapsed)
	fmt.Println("Done")

	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	printDebug(params, ciphertextBefore, valuesWant, ctx.decryptor, ctx.encoder)

	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertextAfter, valuesWant, ctx.decryptor, ctx.encoder)

	return elapsed
}

func sampleRandomInts(gbfv gbfvParameters) ([]*big.Int, error) {
	randomInts := make([]*big.Int, gbfv.slots)
	for i := range randomInts {
		value, err := rand.Int(rand.Reader, gbfv.upperBound)
		if err != nil {
			return nil, err
		}
		randomInts[i] = value
	}
	return randomInts, nil
}

func gbfvCoefficientValues(randomInts []*big.Int, gbfv gbfvParameters) []*big.Float {
	values := make([]*big.Float, gbfv.ringDegree)
	for i := range values {
		idx := int64(i) % gbfv.slots
		block := int64(i) / gbfv.slots
		exponent := big.NewInt(gbfv.ringDegree/gbfv.slots - block - 1)

		coefficient := new(big.Int).Exp(gbfv.baseInt, exponent, nil)
		coefficient.Mul(coefficient, randomInts[idx])
		coefficient.Mod(coefficient, gbfv.upperBound)

		unit := new(big.Float).SetPrec(gbfv.precision).SetInt(coefficient)
		unit.Quo(unit, gbfv.upperBoundFloat)

		values[i] = new(big.Float).SetPrec(gbfv.precision).Mul(unit, gbfv.qFloat)
	}
	return values
}

func normalizeByQ(values []*big.Float, gbfv gbfvParameters) []*big.Float {
	normalized := make([]*big.Float, len(values))
	for i, value := range values {
		normalized[i] = new(big.Float).SetPrec(gbfv.precision).Quo(value, gbfv.qFloat)
	}
	return normalized
}

func perturbValues(values []*big.Float, gbfv gbfvParameters) []*big.Float {
	perturbed := make([]*big.Float, len(values))
	for i, value := range values {
		perturbation := bignum.NewFloat(sampling.RandFloat64(-perturbationMagnitude, perturbationMagnitude), gbfv.precision)
		perturbation.Mul(perturbation, gbfv.qFloat)

		perturbed[i] = new(big.Float).SetPrec(gbfv.precision).Add(value, perturbation)
	}
	return perturbed
}

func inverseTMultiplierValues(gbfv gbfvParameters) []*big.Float {
	values := make([]*big.Float, gbfv.ringDegree)
	for i := range values {
		values[i] = new(big.Float).SetPrec(gbfv.precision)
		if int64(i)%gbfv.slots != 0 {
			continue
		}

		block := int64(i) / gbfv.slots
		exponent := big.NewInt(gbfv.ringDegree/gbfv.slots - block - 1)

		coefficient := new(big.Int).Exp(gbfv.baseInt, exponent, nil)
		unit := new(big.Float).SetPrec(gbfv.precision).SetInt(coefficient)
		unit.Quo(unit, gbfv.upperBoundFloat)

		values[i].Mul(unit, gbfv.qPrimeFloat)
	}
	return values
}

func tPolynomialBigFloat(size int, gbfv gbfvParameters) []*big.Float {
	values := make([]*big.Float, size)
	for i := range values {
		values[i] = new(big.Float).SetPrec(gbfv.precision)
	}
	values[0].SetInt64(gbfv.base)
	values[gbfv.slots].SetInt64(-1)
	return values
}

func encodePlaintexts(params ckks.Parameters, encoder *ckks.Encoder, gbfv gbfvParameters, values []*big.Float) (plaintext, ptT, ptInverseT *rlwe.Plaintext, err error) {
	inputLevel := params.LevelsConsumedPerRescaling() - 1
	inverseTLevel := 2*params.LevelsConsumedPerRescaling() - 1

	plaintext = ckks.NewPlaintext(params, inputLevel)
	plaintext.IsBatched = false
	if err = encoder.EncodePoly(values, plaintext); err != nil {
		return nil, nil, nil, err
	}

	ptT = ckks.NewPlaintext(params, inputLevel)
	ptT.IsBatched = false
	if err = encoder.EncodePoly(tPolynomialBigFloat(int(gbfv.ringDegree), gbfv), ptT); err != nil {
		return nil, nil, nil, err
	}
	ptT.Scale = rlwe.NewScale(gbfv.oneFloat)

	ptInverseT = ckks.NewPlaintext(params, inverseTLevel)
	ptInverseT.IsBatched = false
	if err = encoder.EncodePoly(inverseTMultiplierValues(gbfv), ptInverseT); err != nil {
		return nil, nil, nil, err
	}
	ptInverseT.Scale = rlwe.NewScale(gbfv.qPrime)

	return plaintext, ptT, ptInverseT, nil
}

func bootstrapGBFV(eval *bootstrapping.Evaluator, ciphertext *rlwe.Ciphertext, ptT, ptInverseT *rlwe.Plaintext, gbfv gbfvParameters) (ciphertextBefore, ciphertextAfter *rlwe.Ciphertext, err error) {
	ctTimesT, err := eval.MulNew(ciphertext, ptT)
	if err != nil {
		return nil, nil, err
	}

	ciphertextSub := ciphertext.CopyNew()
	ciphertext.Scale = rlwe.NewScale(gbfv.qFloat)

	fmt.Println("Bootstrapping...")
	bootstrapped, err := eval.Bootstrap(ctTimesT)
	if err != nil {
		return nil, nil, err
	}

	correction, err := eval.MulNew(bootstrapped, ptInverseT)
	if err != nil {
		return nil, nil, err
	}
	if err = eval.Rescale(correction, correction); err != nil {
		return nil, nil, err
	}
	if err = eval.Rescale(correction, correction); err != nil {
		return nil, nil, err
	}

	ciphertextAfter, err = eval.SubNew(ciphertextSub, correction)
	if err != nil {
		return nil, nil, err
	}
	ciphertextAfter.Scale = rlwe.NewScale(gbfv.qFloat)

	return ciphertext, ciphertextAfter, nil
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []*big.Float, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []*big.Float) {
	prec := params.EncodingPrecision()
	valuesTest = make([]*big.Float, params.N())
	valuesWantAdjusted := make([]*big.Float, params.N())
	one := new(big.Float).SetPrec(prec).SetInt64(1)

	for i := range valuesTest {
		valuesTest[i] = new(big.Float).SetPrec(prec)
		valuesWantAdjusted[i] = new(big.Float).SetPrec(prec)
	}

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	for i := range valuesTest {
		valuesTest[i], valuesWantAdjusted[i] = adjustUnitWrap(valuesTest[i], valuesWant[i], one, prec)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.27f %6.27f...\n", valuesTest[0], valuesTest[1])
	fmt.Printf("ValuesWant: %6.27f %6.27f...\n", valuesWantAdjusted[0], valuesWantAdjusted[1])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWantAdjusted, valuesTest, 0, false)
	fmt.Println(precStats.String())
	fmt.Println()

	return valuesTest
}

func adjustUnitWrap(valueTest, valueWant, one *big.Float, precision uint) (adjustedTest, adjustedWant *big.Float) {
	adjustedTest = new(big.Float).SetPrec(precision).Set(valueTest)
	adjustedWant = new(big.Float).SetPrec(precision).Set(valueWant)

	diff := new(big.Float).SetPrec(precision).Sub(adjustedTest, adjustedWant)
	switch {
	case diff.Cmp(big.NewFloat(0.5)) == 1:
		adjustedWant.Add(adjustedWant, one)
	case diff.Cmp(big.NewFloat(-0.5)) == -1:
		adjustedTest.Add(adjustedTest, one)
	}

	return adjustedTest, adjustedWant
}
