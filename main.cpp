#include <cstddef>
#include <sys/time.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <gmp.h>
#include <omp.h>
#include <paillier.h>

#define VEC_SIZE 1000

// Simple class to measure time for each method
class Timer
{
public:
    void start() { m_start = my_clock(); }
    void stop() { m_stop = my_clock(); }
    double elapsed_time() const {
        return m_stop - m_start;
    }

private:
    double m_start, m_stop;
    double my_clock() const {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
};

void noPackingSum(long u[], long v[], long result[], FHESecKey sk, FHEPubKey pk)
{
    Timer tNoPackingEncryption;
    tNoPackingEncryption.start();
    // Vectors to hold the ciphertexts created from the elements of u and v
    std::vector<Ctxt> encU(VEC_SIZE, Ctxt(pk));
    std::vector<Ctxt> encV(VEC_SIZE, Ctxt(pk));

    // Each element is encrypted individually
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	Ctxt tempU(pk);
	pk.Encrypt(tempU, to_ZZX(u[i]));
	encU[i] = tempU;

	Ctxt tempV(pk);
	pk.Encrypt(tempV, to_ZZX(v[i]));
	encV[i] = tempV;
    }
    tNoPackingEncryption.stop();
    std::cout << "HElib No packing encryption: " << tNoPackingEncryption.elapsed_time() << "s." <<  std::endl;

    Timer tNoPackingSum;
    tNoPackingSum.start();
    // Sum vectors element wise
    // The sum is stored in U
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	encU[i] += encV[i];
    }
    tNoPackingSum.stop();
    std::cout << "HElib No packing sum: " << tNoPackingSum.elapsed_time() << "s." <<  std::endl;

    Timer tNoPackingDecryption;
    tNoPackingDecryption.start();
    // Decrypt result
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	ZZX element;
	sk.Decrypt(element, encU[i]);
	result[i] = conv<long>(element[0]);
    }
    tNoPackingDecryption.stop();
    std::cout << "HElib No packing decryption: " << tNoPackingDecryption.elapsed_time() << "s." <<  std::endl;
}

void packingPolySum(long u[], long v[], long result[], FHESecKey sk, FHEPubKey pk)
{
    Timer tPackingPolyEncryption;
    tPackingPolyEncryption.start();
    // ZZX is a class for polynomials from the NTL library
    ZZX U, V;                               

    // Set the length of the polynomial U(x) and V(x)
    U.SetLength(VEC_SIZE);
    V.SetLength(VEC_SIZE);

    // Set the coefficients of the polynomials U(x) and V(x).
    for (int i = 0; i < VEC_SIZE; i++) {
	SetCoeff(U, i, u[i]);
	SetCoeff(V, i, v[i]);	
    }

    // Ciphertexts that will hold the polynomials encrypted using public key pk
    Ctxt encU(pk);                          
    Ctxt encV(pk);                          

    // Encrypt the polynomials into the ciphertexts
    pk.Encrypt(encU, U);
    pk.Encrypt(encV, V);
    tPackingPolyEncryption.stop();
    std::cout << "HElib poly packing encryption: " << tPackingPolyEncryption.elapsed_time() << "s." <<  std::endl;

    Timer tPackingPolySum;
    tPackingPolySum.start();
    // Multiply the ciphertexts and store the result into encU
    encU += encV;
    tPackingPolySum.stop();
    std::cout << "HElib poly packing sum: " << tPackingPolySum.elapsed_time() << "s." <<  std::endl;    

    Timer tPackingPolyDecryption;
    tPackingPolyDecryption.start();
    // Decrypt the multiplied ciphertext into a polynomial using the secret key sk
    ZZX resultPoly;
    sk.Decrypt(resultPoly, encU);

    // Assign the values of the polynomial's coefficients to the result vector
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	result[i] = conv<long>(resultPoly[i]);
    }
    tPackingPolyDecryption.stop();
    std::cout << "HElib poly packing decryption: " << tPackingPolyDecryption.elapsed_time() << "s." <<  std::endl;
}

void packingSubfieldSum(long u[], long v[], long result[], FHEPubKey pk, FHESecKey sk, FHEcontext& context)
{
    Timer tPackingSubfieldEncryption;
    tPackingSubfieldEncryption.start();
    // Creates a helper object based on the context
    EncryptedArray ea(context, context.alMod.getFactorsOverZZ()[0]); 

    // The vectors should have the same size as the EncryptedArray (ea.size),
    // so fill the other positions with 0 which won't change the result
    std::vector<long> U(u, u + VEC_SIZE);
    std::vector<long> V(v, v + VEC_SIZE);
    for (int i = VEC_SIZE; i < ea.size(); i++) {
	U.push_back(0);
	V.push_back(0);
    }

    // Ciphertexts that will hold the encrypted vectors
    Ctxt encU(pk);
    Ctxt encV(pk);

    // Encrypt the whole vector into one ciphertext using packing
    ea.encrypt(encU, pk, U);
    ea.encrypt(encV, pk, V);
    tPackingSubfieldEncryption.stop();
    std::cout << "HElib subfield packing encryption: " << tPackingSubfieldEncryption.elapsed_time() << "s." <<  std::endl;    

    Timer tPackingSubfieldSum;
    tPackingSubfieldSum.start();
    // Sum ciphertexts and set the result to encU
    encU += encV;
    tPackingSubfieldSum.stop();
    std::cout << "HElib subfield packing sum: " << tPackingSubfieldSum.elapsed_time() << "s." <<  std::endl;        

    Timer tPackingSubfieldDecryption;
    tPackingSubfieldDecryption.start();
    // Decrypt the result
    std::vector<long> res(ea.size(), 0);
    ea.decrypt(encU, sk, res);   

    // Assign the values of the polynomial's coefficients to the result vector
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	result[i] = conv<long>(res[i]);
    }
    tPackingSubfieldDecryption.stop();
    std::cout << "HElib subfield packing decryption: " << tPackingSubfieldDecryption.elapsed_time() << "s." <<  std::endl;            
}

void paillierSum(long u[], long v[], long result[], paillier_pubkey_t* pubKey, paillier_prvkey_t* secKey)
{
    Timer tPaillierEncryption;
    tPaillierEncryption.start();
    // Vectors that will hold the ciphertexts
    std::vector<paillier_ciphertext_t*> encU(VEC_SIZE);
    std::vector<paillier_ciphertext_t*> encV(VEC_SIZE);    

    // Encrypt elements
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	paillier_plaintext_t* ptxtU = paillier_plaintext_from_ui((int)u[i]);
	paillier_ciphertext_t* ctxtU = paillier_enc(NULL, pubKey, ptxtU, paillier_get_rand_devurandom);
	encU[i] = ctxtU;


	paillier_plaintext_t* ptxtV = paillier_plaintext_from_ui((int)v[i]);
	paillier_ciphertext_t* ctxtV = paillier_enc(NULL, pubKey, ptxtV, paillier_get_rand_devurandom);
	encV[i] = ctxtV;

	paillier_freeplaintext(ptxtU);		
	paillier_freeplaintext(ptxtV);	
    }
    tPaillierEncryption.stop();
    std::cout << "Paillier encryption: " << tPaillierEncryption.elapsed_time() << "s." <<  std::endl;                

    Timer tPaillierSum;
    tPaillierSum.start();
    // Sum encrypted vectors element wise
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	paillier_ciphertext_t* encryptedSum = paillier_create_enc_zero();
	paillier_mul(pubKey, encryptedSum, encU[i], encV[i]);
	paillier_freeciphertext(encU[i]);
	encU[i] = encryptedSum;
    }
    tPaillierSum.stop();
    std::cout << "Paillier sum: " << tPaillierSum.elapsed_time() << "s." <<  std::endl;                

    Timer tPaillierDecryption;
    tPaillierDecryption.start();
    // Decrypt the result
    #pragma omp parallel for
    for (int i = 0; i < VEC_SIZE; i++) {
	paillier_plaintext_t* dec;
	dec = paillier_dec(NULL, pubKey, secKey, encU[i]);
	result[i] = mpz_get_ui(dec->m);
	
	paillier_freeciphertext(encU[i]);
	paillier_freeciphertext(encV[i]);
	paillier_freeplaintext(dec);	
    }
    tPaillierDecryption.stop();
    std::cout << "Paillier decryption: " << tPaillierDecryption.elapsed_time() << "s." <<  std::endl;                
}

int main(int argc, char **argv)
{
    /*** BEGIN INITIALIZATION ***/
    long m = 0;                   // Specific modulus
    long p = 257;                 // Plaintext base [default=2], should be a prime number
    long r = 3;                   // Lifting [default=1]
    long L = 10;                  // Number of levels in the modulus chain [default=heuristic]
    long c = 2;                   // Number of columns in key-switching matrix [default=2]
    long w = 64;                  // Hamming weight of secret key
    long d = 1;                   // Degree of the field extension [default=1]
    long k = 80;                  // Security parameter [default=80] 
    long s = 0;                   // Minimum number of slots [default=0]
    
    Timer tInit;
    tInit.start();
	
    std::cout << "Finding m... " << std::flush;
    m = FindM(k, L, c, p, d, s, 0);           // Find a value for m given the specified values
    std::cout << "m = " << m << std::endl;
	
    std::cout << "Initializing context... " << std::flush;
    FHEcontext context(m, p, r); 	          // Initialize context
    buildModChain(context, L, c);             // Modify the context, adding primes to the modulus chain
    std::cout << "OK!" << std::endl;

    std::cout << "Generating keys... " << std::flush;
    FHESecKey sk(context);                    // Construct a secret key structure
    const FHEPubKey& pk = sk;                 // An "upcast": FHESecKey is a subclass of FHEPubKey
    sk.GenSecKey(w);                          // Actually generate a secret key with Hamming weight
    addSome1DMatrices(sk);                    // Extra information for relinearization
    std::cout << "OK!" << std::endl;

    // Arrays whose elements will be the coefficients of the polynomials U(x) and V(x)
    long u[VEC_SIZE];
    long v[VEC_SIZE];
    long result[VEC_SIZE];    

    // Initialize arrays
    // Array u will have even numbers, and array v will have odd numbers
    for (int i = 0; i < VEC_SIZE; i++) {
	u[i] = 2*i;
	v[i] = 2*i + 1;
    }
    tInit.stop();
    std::cout << "Time taken for the initialization: " << tInit.elapsed_time() << std::endl;
	
    /*** METHOD 1: SUM ARRAYS USING HELIB WITHOUT PACKING ***/
    Timer tMethod1;
    tMethod1.start();
    noPackingSum(u, v, result, sk, pk);
    tMethod1.stop();
    std::cout << "HElib without packing method done in " << tMethod1.elapsed_time() << "s." <<  std::endl;	
   
    /*** METHOD 2: SUM ARRAYS USING HELIB WITH POLYNOMIAL PACKING***/
    Timer tMethod2;
    tMethod2.start();
    packingPolySum(u, v, result, sk, pk);
    tMethod2.stop();
    std::cout << "HElib with polynomial packing method done in " << tMethod2.elapsed_time() << "s." <<  std::endl;

   
    /*** METHOD 3: SUM ARRAYS USING HELIB WITH SUBFIELD PACKING ***/
    Timer tMethod3;
    tMethod3.start();
    packingSubfieldSum(u, v, result, pk, sk, context);
    tMethod3.stop();
    std::cout << "HElib with subfield packing method done in " << tMethod3.elapsed_time() << "s." <<  std::endl;

    /*** METHOD 4: SUM ARRAYS USING PAILLIER LIB***/
    Timer tMethod4;
    tMethod4.start();

    Timer tPaillierInitialization;
    tPaillierInitialization.start();
    // Security parameter (number of bits of the modulus)
    const long n = 256;   
    
    // Generate public and secret keys
    paillier_pubkey_t* pubKey;
    paillier_prvkey_t* secKey;
    paillier_keygen(n, &pubKey, &secKey, paillier_get_rand_devurandom);
    tPaillierInitialization.stop();
    std::cout << "Paillier initialization: " << tPaillierInitialization.elapsed_time() << "s." <<  std::endl;    
    
    paillierSum(u, v, result, pubKey, secKey);

    paillier_freepubkey(pubKey);
    paillier_freeprvkey(secKey);
    
    tMethod4.stop();
    std::cout << "Paillier lib method done in " << tMethod4.elapsed_time() << "s." <<  std::endl;

    // // Checking results
    // for (int i = 0; i < VEC_SIZE; i++) {
    // 	std::cout << u[i] << "+" << v[i] << " = " << result[i] << std::endl;
    // }

    return 0;
}
