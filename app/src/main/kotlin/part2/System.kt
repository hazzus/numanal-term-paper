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
    K1,
    K2,
    K3,
    D_HCl,
    D_H2,
    D_AlCl2,
    D_AlCl3,
    D_AlCl,
    P_AlCl,
    P_H2,
    P_HCl,
    P_AlCl2,
    P_AlCl3,
    PG_AlCl,
    PG_HCl,
    PG_H2,
    PG_AlCl2,
    PG_AlCl3
}

val part2System: Env<Part2Coefficient>.(Vector) -> Vector = intercept {
    listOf(
        K1 * P_AlCl.pow(2.0) * P_H2 - P_HCl.pow(2.0),
        K2 * P_AlCl2 * P_H2 - P_HCl.pow(2.0),
        K3 * P_AlCl3.pow(2.0) * P_H2.pow(3.0) - P_HCl.pow(6),
        D_HCl * (PG_HCl - P_HCl) + 2.0 * D_H2 * (PG_H2 - P_H2),
        D_AlCl * (PG_AlCl - P_AlCl) + 2.0 * D_AlCl2 * (PG_AlCl2 - P_AlCl2) + 3.0 * D_AlCl3 * (PG_AlCl3 - P_AlCl3) + D_HCl * (PG_HCl - P_HCl)
    )
}

val part2Jacobi: Matrix<Env<Part2Coefficient>.(Vector) -> Double> = listOf(
    interceptedListOf(
        { 2.0 * K1 * P_AlCl * P_H2 }, // alcl
        { K1 * P_AlCl.pow(2.0) }, // h2
        { -2.0 * P_HCl }, // hcl
        zero, // alcl2
        zero // alcl3
    ),
    interceptedListOf(
        zero,
        { K2 * P_AlCl2 },
        { -2.0 * P_HCl },
        { K2 * P_H2 },
        zero
    ),
    interceptedListOf(
        zero,
        { 3.0 * K3 * P_AlCl3.pow(2.0) * P_H2.pow(2.0) },
        { -6.0 * P_HCl.pow(5.0) },
        zero,
        { K3 * 2.0 * P_AlCl3 * P_H2.pow(3.0) }),
    interceptedListOf(
        zero,
        { -2.0 * D_H2 },
        { D_HCl * (-1.0) },
        zero,
        zero
    ),
    interceptedListOf(
        { D_AlCl * (-1.0) },
        zero,
        { (-1.0) * D_HCl },
        { -2.0 * D_AlCl2 },
        { -3.0 * D_AlCl3 }
    )
)

class Part2Env(private val coeffs: Coefficients) : Env<Part2Coefficient> {
    override val Part2Coefficient.c: Double
        get() = when (this) {
            K1 -> coeffs.K4
            K2 -> coeffs.K5
            K3 -> coeffs.K6
            D_HCl -> coeffs.dHCl
            D_H2 -> coeffs.dH2
            D_AlCl -> coeffs.dGaCl
            D_AlCl2 -> coeffs.dGaCl2
            D_AlCl3 -> coeffs.dGaCl3
            PG_H2 -> coeffs.pGH2
            PG_HCl -> coeffs.pGHCl
            PG_AlCl -> coeffs.pGGaCl
            PG_AlCl2 -> coeffs.pGGaCl2
            PG_AlCl3 -> coeffs.pGGaCl3
            else -> error(":(")
        }

    override val variableIdentifiers: List<Part2Coefficient> = listOf(P_AlCl, P_H2, P_HCl, P_AlCl2, P_AlCl3)
}

