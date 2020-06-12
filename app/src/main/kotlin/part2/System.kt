package part2

import Env
import Matrix
import Vector
import intercept
import interceptedListOf
import part2.Part2Coefficient.*
import zero

@Suppress("EnumEntryName")
enum class Part2Coefficient {
    K4,
    K5,
    K6,
    D_HCl,
    D_H2,
    D_GaCl2,
    D_GaCl3,
    D_GaCl,
    P_GaCl,
    P_H2,
    P_HCl,
    P_GaCl2,
    P_GaCl3,
    PG_GaCl,
    PG_HCl,
    PG_H2,
    PG_GaCl2,
    PG_GaCl3
}

val part2System: Env<Part2Coefficient>.(Vector) -> Vector = intercept {
    listOf(
        K4 * P_GaCl.pow(2.0) * P_H2 - P_HCl.pow(2.0),
        K5 * P_GaCl2 * P_H2 - P_HCl.pow(2.0),
        K6 * P_GaCl3.pow(2.0) * P_H2.pow(3.0) - P_HCl.pow(6),
        D_HCl * (PG_HCl - P_HCl) + 2.0 * D_H2 * (PG_H2 - P_H2),
        D_GaCl * (PG_GaCl - P_GaCl) + 2.0 * D_GaCl2 * (PG_GaCl2 - P_GaCl2) + 3.0 * D_GaCl3 * (PG_GaCl3 - P_GaCl3) + D_HCl * (PG_HCl - P_HCl)
    )
}

val part2Jacobi: Matrix<Env<Part2Coefficient>.(Vector) -> Double> = listOf(
    interceptedListOf(
        { 2.0 * K4 * P_GaCl * P_H2 }, // gacl
        { K5 * P_GaCl.pow(2.0) }, // h2
        { -2.0 * P_HCl }, // hcl
        zero, // gacl2
        zero // gacl3
    ),
    interceptedListOf(
        zero,
        { K5 * P_GaCl2 },
        { -2.0 * P_HCl },
        { K5 * P_H2 },
        zero
    ),
    interceptedListOf(
        zero,
        { 3.0 * K6 * P_GaCl3.pow(2.0) * P_H2.pow(2.0) },
        { -6.0 * P_HCl.pow(5.0) },
        zero,
        { K6 * 2.0 * P_GaCl3 * P_H2.pow(3.0) }),
    interceptedListOf(
        zero,
        { -2.0 * D_H2 },
        { D_HCl * (-1.0) },
        zero,
        zero
    ),
    interceptedListOf(
        { D_GaCl * (-1.0) },
        zero,
        { (-1.0) * D_HCl },
        { -2.0 * D_GaCl2 },
        { -3.0 * D_GaCl3 }
    )
)

class Part2Env(private val coeffs: Coefficients) : Env<Part2Coefficient> {
    override val Part2Coefficient.c: Double
        get() = when (this) {
            K4 -> coeffs.K4
            K5 -> coeffs.K5
            K6 -> coeffs.K6
            D_HCl -> coeffs.dHCl
            D_H2 -> coeffs.dH2
            D_GaCl -> coeffs.dGaCl
            D_GaCl2 -> coeffs.dGaCl2
            D_GaCl3 -> coeffs.dGaCl3
            PG_H2 -> coeffs.pGH2
            PG_HCl -> coeffs.pGHCl
            PG_GaCl -> coeffs.pGGaCl
            PG_GaCl2 -> coeffs.pGGaCl2
            PG_GaCl3 -> coeffs.pGGaCl3
            else -> error(":(")
        }

    override val variableIdentifiers: List<Part2Coefficient> = listOf(P_GaCl, P_H2, P_HCl, P_GaCl2, P_GaCl3)
}

