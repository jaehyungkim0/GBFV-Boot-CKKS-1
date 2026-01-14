// GBFV Bootstrapping from CKKS, relying on META-BTS to instantiate high precision CKKS bootstrapping

// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"crypto/rand"
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

func main() {

	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16

	if *flagShort {
		LogN -= 3
	}

	//==============================
	//=== 1) RESIDUAL PARAMETERS ===
	//==============================

	// First we must define the residual parameters.
	// The residual parameters are the parameters used outside of the bootstrapping circuit.
	// For this example, we have a LogN=16, logQ = (55+45) + 5*(45+45) and logP = 3*61, so LogQP = 638.
	// With LogN=16, LogQP=638 and H=192, these parameters achieve well over 128-bit of security.
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,              // Log2 of the ring degree
		LogQ:            []int{60, 45, 45, 45},     // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61, 61}, // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 90,                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	prec := params.EncodingPrecision()

	//==========================================
	//=== 2) BOOTSTRAPPING PARAMETERSLITERAL ===
	//==========================================

	// The bootstrapping circuit use its own Parameters which will be automatically
	// instantiated given the residual parameters and the bootstrapping parameters.

	// !WARNING! The bootstrapping parameters are not ensure to be 128-bit secure, it is the
	// responsibility of the user to check that the meet the security requirement and tweak them if necessary.

	// Note that the default bootstrapping parameters use LogN=16 and a ternary secret with H=192 non-zero coefficients
	// which provides parameters which are at least 128-bit if their LogQP <= 1550.

	// For this first example, we do not specify any circuit specific optional field in the bootstrapping parameters literal.
	// Thus we expect the bootstrapping to give an average precision of 27.9 bits with H=192 (and 24.4 with H=N/2)
	// if the plaintext values are uniformly distributed in [-1, 1] for both the real and imaginary part.
	// See `circuits/bootstrapping/parameters_literal.go` for detailed information about the optional fields.
	btpParametersLit := bootstrapping.ParametersLiteral{
		// We specify LogN to ensure that both the residual parameters and the bootstrapping parameters
		// have the same LogN. This is not required, but we want it for this example.
		LogN: utils.Pointy(LogN),

		// In this example we need manually specify the number of auxiliary primes (i.e. #Pi) used by the
		// evaluation keys of the bootstrapping circuit, so that the size of LogQP  meets the security target.
		LogP: []int{61, 61, 61, 61},

		// Sets the IterationsParameters.
		// The default bootstrapping parameters have 27.9 bits of average precision and
		// ~25 bits of minimum precision, and the maximum precision that can be theoretically
		// achieved is LogScale - LogN/2.
		// Therefore we start with 27.9 bits and each can in theory increase the precision an additional 25 bits.
		// However, to achieve the best possible precision, we must carefully adjust each iteration by hand so
		// that the sum of all the minimum precision is as close as possible
		// to LogScale - LogN/2. Here 27.9+25+25+5 ~= 82.5 (for the insecure parameters with LogN=13, with
		// the secure parameters using LogN=16 achieve 82.5 - (16-13)/2 = 81 bits of precision).
		IterationsParameters: &bootstrapping.IterationsParameters{
			BootstrappingPrecision: []float64{25, 25, 5},
			ReservedPrimeBitSize:   28,
		},

		// In this example we manually specify the bootstrapping parameters' secret distribution.
		// This is not necessary, but we ensure here that they are the same as the residual parameters.
		Xs: params.Xs(),
	}

	//===================================
	//=== 3) BOOTSTRAPPING PARAMETERS ===
	//===================================

	// Now that the residual parameters and the bootstrapping parameters literals are defined, we can instantiate
	// the bootstrapping parameters.
	// The instantiated bootstrapping parameters store their own ckks.Parameter, which are the parameters of the
	// ring used by the bootstrapping circuit.
	// The bootstrapping parameters are a wrapper of ckks.Parameters, with additional information.
	// They therefore has the same API as the ckks.Parameters and we can use this API to print some information.
	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if *flagShort {
		// Corrects the message ratio Q0/|m(X)| to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	// We print some information about the residual parameters.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(), params.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	// And some information about the bootstrapping parameters.
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//===========================
	//=== 4) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//========================
	//=== 5) BOOTSTRAPPING ===
	//========================

	// Instantiates the bootstrapper
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}
	
	// Paramters
	// We are using t(X) = X^k - b over Z[X]/(X^N+1)
	ringQ := params.RingQ().AtLevel(params.LevelsConsumedPerRescaling()-1)
	Q := ringQ.ModulusAtLevel[params.LevelsConsumedPerRescaling()-1]
	Qprime := new(big.Int).Div(ringQ.ModulusAtLevel[2 * params.LevelsConsumedPerRescaling() - 1], Q)
	fmt.Println("Q = ", Q)
	f := new(big.Float).SetPrec(prec)
	f.SetInt(Q)
	fprime := new(big.Float).SetPrec(prec)
	fprime.SetInt(Qprime)
	b := int64(256)
	bbig := big.NewInt(b)
	k := int64(1)
	N := (int64)(1 << LogN)
	one := new(big.Float).SetPrec(prec)
	one.SetInt(big.NewInt(1))
	// Compute b^{N/k}
	upperBound := new(big.Int).Exp(bbig, big.NewInt(N/k), nil)

        // Add 1: upperBound = b^{N/k} + 1
        upperBound.Add(upperBound, big.NewInt(1))
	upperBoundFloat := new(big.Float).SetPrec(prec)
	upperBoundFloat.SetInt(upperBound)

	fmt.Printf("GBFV Bootstrapping Parameters: b = %d, k = %d\n", b, k) 


	var totalDuration time.Duration
	num_iter := 20

	for iter := 0; iter < num_iter; iter++ {
	fmt.Println("-----------------------------------------------------")
	fmt.Printf("------------------ Iteration %d ----------------------\n", iter)
	fmt.Println("-----------------------------------------------------")

        // Generate k random integer in [0, b^{N/k} + 1)
	randInts := make([]*big.Int, k)
	for i := range randInts {
        	randInts[i], err = rand.Int(rand.Reader, upperBound)
        	if err != nil {
        		panic(err)
        	}
	}

	valuesWant := make([]*big.Float, N)

	// Generate a corresponding GBFV ciphertext
	for i := range valuesWant {
		valuesWant[i] = new(big.Float).SetPrec(prec)
		idx := int64(i) % k
		j := int64(i) / k
		Nmj := big.NewInt(N/k-int64(j)-1)
		bint := new(big.Int).Exp(bbig, Nmj, nil)
		bint = new(big.Int).Mul(bint, randInts[idx]) // multiply by random integer in [0, b^{N/k}+1)
		bint = new(big.Int).Mod(bint, upperBound) // modular reduce by b^{N/k}+1
		bi := new(big.Float).SetPrec(prec)
		bi.SetInt(bint)
		bi.Quo(bi, upperBoundFloat)
		valuesWant[i].Mul(bi, f)
	}

	valuesWantPert := make([]*big.Float, N)

	// Providing perturbations to generate error	
	for i := range valuesWantPert {
		pert := bignum.NewFloat(sampling.RandFloat64(-0.000000001, 0.000000001), prec)
		pert.Mul(pert, f)
		valuesWantPert[i] = new(big.Float).SetPrec(prec)
		valuesWantPert[i].Add(valuesWant[i], pert)
	}

	for i := range valuesWant {
		valuesWant[i].Quo(valuesWant[i], f)
	}

	// Polynomial -q'/t(X) 
	valuesMult := make([]*big.Float, N)
	for i := range valuesMult {
		valuesMult[i] = new(big.Float).SetPrec(prec)
		if (int64(i) % k == 0) {
			j := int64(i) / k
			Nmj := big.NewInt(N/k-int64(j)-1)
			bint := new(big.Int).Exp(bbig, Nmj, nil)
			bi := new(big.Float).SetPrec(prec)
			bi.SetInt(bint)
			bi.Quo(bi, upperBoundFloat)
			valuesMult[i].Mul(bi, fprime)
		} else {
			valuesMult[i].SetInt(big.NewInt(0))
		}
	}

	// Polynomial -t(X) = b - X^k
	vecpoly := make([]float64, N)
	for i := range vecpoly {
		if i == 0 {
			vecpoly[i] = float64(b)
		} else if int64(i) == k {
			vecpoly[i] = -1.0
		} else {
			vecpoly[i] = 0.0
		}
	}

	// We encrypt at level=LevelsConsumedPerRescaling-1
	plaintext := ckks.NewPlaintext(params, params.LevelsConsumedPerRescaling()-1)
	plaintext.IsBatched = false
	if err := encoder.EncodePoly(valuesWantPert, plaintext); err != nil {
		panic(err)
	}

	pt_mul := ckks.NewPlaintext(params, params.LevelsConsumedPerRescaling()-1)
	pt_mul.IsBatched = false
	if err := encoder.EncodePoly(vecpoly, pt_mul); err != nil {
		panic(err)
	}
	pt_mul.Scale = rlwe.NewScale(one);

	plaintextMult := ckks.NewPlaintext(params, 2 * params.LevelsConsumedPerRescaling()-1)
	plaintextMult.IsBatched = false
	if err := encoder.EncodePoly(valuesMult, plaintextMult); err != nil {
		panic(err)
	}
	plaintextMult.Scale = rlwe.NewScale(Qprime)

	// Encrypt
	ciphertext1, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	start := time.Now()
	// Multiply -t(X)
	ciphertext3, err := eval.MulNew(ciphertext1, pt_mul)
	ciphertextSub := ciphertext1.CopyNew()

	ciphertext1.Scale = rlwe.NewScale(f)

	// CKKS-Bootstrap the ciphertext 
	fmt.Println("Bootstrapping...")
	ciphertext2, err := eval.Bootstrap(ciphertext3)

	// Multiply -q'/t(X)
	ciphertext2, err = eval.MulNew(ciphertext2, plaintextMult)
	err = eval.Rescale(ciphertext2, ciphertext2)
	err = eval.Rescale(ciphertext2, ciphertext2)

	// Subtraction Step	
	ciphertext2, err = eval.SubNew(ciphertextSub, ciphertext2)

	elapsed := time.Since(start) // measure elapsed time
	totalDuration += elapsed
	fmt.Printf("GBFV Bootstrapping Time (iteration %d): %s\n", iter, elapsed)

	ciphertext2.Scale = rlwe.NewScale(f)

	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//==================
	//=== 6) DECRYPT ===
	//==================

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	printDebug(params, ciphertext1, valuesWant, decryptor, encoder)
	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, ciphertext2, valuesWant, decryptor, encoder)

	} // end iteration for loop
	fmt.Println("-------------------------------------------------------------------")
	average := totalDuration / time.Duration(num_iter)
	fmt.Println("Average Bootstrapping Time:", average)
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []*big.Float, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []*big.Float) {
	prec := params.EncodingPrecision()
	valuesTest = make([]*big.Float, 1 << params.LogN())
	valuesWant2 := make([]*big.Float, 1 << params.LogN())
	one := new(big.Float).SetPrec(prec)
	one.SetInt(big.NewInt(1))
	for i := range valuesTest {
		valuesTest[i] = new(big.Float).SetPrec(prec)
		valuesWant2[i] = new(big.Float).SetPrec(prec)
	}

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	for i := range valuesTest {
		temp := new(big.Float).SetPrec(prec)
		temp.Sub(valuesTest[i], valuesWant[i])
		
		if temp.Cmp(big.NewFloat(0.5)) == 1 {
			valuesWant2[i].Add(valuesWant[i], one)
		} else if temp.Cmp(big.NewFloat(-0.5)) == -1 {
			valuesTest[i].Add(valuesTest[i], one)
			valuesWant2[i] = valuesWant[i]
		} else {
			valuesWant2[i] = valuesWant[i]
		}
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.27f %6.27f...\n", valuesTest[0], valuesTest[1])
	fmt.Printf("ValuesWant: %6.27f %6.27f...\n", valuesWant[0], valuesWant[1])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant2, valuesTest, 0, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
