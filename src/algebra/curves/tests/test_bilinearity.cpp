/**
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <iostream>
#include "common/profiling.hpp"
//#include "algebra/curves/edwards/edwards_pp.hpp"
#ifdef CURVE_BN128
#include "algebra/curves/bn128/bn128_pp.hpp"
#endif
#include "algebra/curves/alt_bn128/alt_bn128_pp.hpp"
//#include "algebra/curves/mnt/mnt4/mnt4_pp.hpp"
//#include "algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "algebra/curves/alt_bn128/alt_bn128_pairing.hpp"
#include "algebra/curves/alt_bn128/alt_bn128_pairing.cpp"

using namespace libsnark;

template<typename ppT>
void pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    //G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    //printf("%llu\n", Fr<ppT>::random_element());
    G1<ppT> P = Fr<ppT>("2") * G1<ppT>::one();
    //G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();
    G2<ppT> Q = Fr<ppT>("3") * G2<ppT>::one();

    printf("P:\n");
    P.print();
    P.print_coordinates();
    printf("Q:\n");
    Q.print();
    Q.print_coordinates();
    printf("\n\n");

    //Fr<ppT> s = Fr<ppT>::random_element();
    Fr<ppT> s = Fr<ppT>("2");
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    sP.print_coordinates();
    sQ.print_coordinates();


    alt_bn128_Fq12 elt;
    alt_bn128_ate_G1_precomp y = alt_bn128_ate_precompute_G1(sP);
    alt_bn128_ate_G2_precomp z = alt_bn128_ate_precompute_G2(sQ);
    alt_bn128_Fq12 xx = alt_bn128_ate_miller_loop(y, z);
    std::cout << "ate_miller_loop: " << xx << std::endl;

    //GT<ppT> zz = alt_bn128_final_exponentiation(xx);
    //alt_bn128_Fq12 zzi = zz.inverse();
    //std::cout << "final_exponentiation: " << zz << std::endl;
    GT<ppT> xxi = xx.inverse();
    std::cout << "elt inverse: " << xxi << std::endl;

    alt_bn128_Fq12 A = alt_bn128_Fq12(xx.c0,-xx.c1);
    std::cout << "A: " << A << std::endl;

    alt_bn128_Fq12 C = A * xxi;
    alt_bn128_Fq12 D = C.Frobenius_map(2);
    std::cout << "C: " << C << std::endl;
    std::cout << "D: " << D << std::endl;

    alt_bn128_Fq12 result = D * C;
    std::cout << "result: " << D << std::endl;

    alt_bn128_Fq12 E =  D.cyclotomic_exp(bigint<100>("4965661367192848881"));
    std::cout << "E: " << E << std::endl;

    std::cout << "E unitary_inverse: " << E.unitary_inverse() << std::endl;

    alt_bn128_Fq12 F = alt_bn128_exp_by_neg_z(E);
    std::cout << "F: " << E << std::endl;

    //exit(1);

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::reduced_pairing(P, Q)^s;

    //std::cout << "ans1: " << ans1 << std::endl;

    ans1.print();
    printf("\n\n");
    ans2.print();
    printf("\n\n");
    ans3.print();
    printf("\n\n");
    assert(ans1 == ans2);
    assert(ans2 == ans3);

    assert(ans1 != GT_one);
    assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("\n\n");
}

template<typename ppT>
void double_miller_loop_test()
{
    const G1<ppT> P1 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G1<ppT> P2 = (Fr<ppT>::random_element()) * G1<ppT>::one();
    const G2<ppT> Q1 = (Fr<ppT>::random_element()) * G2<ppT>::one();
    const G2<ppT> Q2 = (Fr<ppT>::random_element()) * G2<ppT>::one();

    const G1_precomp<ppT> prec_P1 = ppT::precompute_G1(P1);
    const G1_precomp<ppT> prec_P2 = ppT::precompute_G1(P2);
    const G2_precomp<ppT> prec_Q1 = ppT::precompute_G2(Q1);
    const G2_precomp<ppT> prec_Q2 = ppT::precompute_G2(Q2);

    const Fqk<ppT> ans_1 = ppT::miller_loop(prec_P1, prec_Q1);
    const Fqk<ppT> ans_2 = ppT::miller_loop(prec_P2, prec_Q2);
    const Fqk<ppT> ans_12 = ppT::double_miller_loop(prec_P1, prec_Q1, prec_P2, prec_Q2);
    assert(ans_1 * ans_2 == ans_12);
}

template<typename ppT>
void affine_pairing_test()
{
    GT<ppT> GT_one = GT<ppT>::one();

    printf("Running bilinearity tests:\n");
    G1<ppT> P = (Fr<ppT>::random_element()) * G1<ppT>::one();
    G2<ppT> Q = (Fr<ppT>::random_element()) * G2<ppT>::one();

    printf("P:\n");
    P.print();
    printf("Q:\n");
    Q.print();
    printf("\n\n");

    Fr<ppT> s = Fr<ppT>::random_element();
    G1<ppT> sP = s * P;
    G2<ppT> sQ = s * Q;

    alt_bn128_Fq12 elt;
    alt_bn128_pp::init_public_params();
    std::cout << "precompute_G1: " << alt_bn128_pp::precompute_G1(P) << std::endl; 

    printf("Pairing bilinearity tests (three must match):\n");
    GT<ppT> ans1 = ppT::affine_reduced_pairing(sP, Q);
    GT<ppT> ans2 = ppT::affine_reduced_pairing(P, sQ);
    GT<ppT> ans3 = ppT::affine_reduced_pairing(P, Q)^s;
    ans1.print();
    ans2.print();
    ans3.print();
    assert(ans1 == ans2);
    assert(ans2 == ans3);

    assert(ans1 != GT_one);
    assert((ans1^Fr<ppT>::field_char()) == GT_one);
    printf("\n\n");
}

int main(void)
{
    start_profiling();
/*
    edwards_pp::init_public_params();
    pairing_test<edwards_pp>();
    double_miller_loop_test<edwards_pp>();

    mnt6_pp::init_public_params();
    pairing_test<mnt6_pp>();
    double_miller_loop_test<mnt6_pp>();
    affine_pairing_test<mnt6_pp>();

    mnt4_pp::init_public_params();
    pairing_test<mnt4_pp>();
    double_miller_loop_test<mnt4_pp>();
    affine_pairing_test<mnt4_pp>();
*/
    alt_bn128_pp::init_public_params();
    pairing_test<alt_bn128_pp>();
    double_miller_loop_test<alt_bn128_pp>();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    bn128_pp::init_public_params();
    pairing_test<bn128_pp>();
    double_miller_loop_test<bn128_pp>();
#endif
}
