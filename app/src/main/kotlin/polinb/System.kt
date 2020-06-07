package polinb

import Env
import Matrix
import Vector
import intercept
import interceptedListOf
import polinb.PolinbCoefficient.*
import zero

@Suppress("EnumEntryName")
enum class PolinbCoefficient {
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

val polinbSystem: Env<PolinbCoefficient>.(Vector) -> Vector = intercept {
    listOf(
        K1 * P_AlCl.pow(2.0) * P_H2 - P_HCl.pow(2.0),
        K2 * P_AlCl2 * P_H2 - P_HCl.pow(2.0),
        K3 * P_AlCl3.pow(2.0) * P_H2.pow(3.0) - P_HCl.pow(6),
        D_HCl * (PG_HCl - P_HCl) + 2.0 * D_H2 * (PG_H2 - P_H2),
        D_AlCl * (PG_AlCl - P_AlCl) + 2.0 * D_AlCl2 * (PG_AlCl2 - P_AlCl2) + 3.0 * D_AlCl3 * (PG_AlCl3 - P_AlCl3) + D_HCl * (PG_HCl - P_HCl)
    )
}

val polinbJacobi: Matrix<Env<PolinbCoefficient>.(Vector) -> Double> = listOf(
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

class PolinbEnv(private val coeffs: Coefficients) : Env<PolinbCoefficient> {
    override val PolinbCoefficient.c: Double
        get() = when (this) {
            K1 -> coeffs.K1
            K2 -> coeffs.K2
            K3 -> coeffs.K3
            D_HCl -> coeffs.dHCl
            D_H2 -> coeffs.dH2
            D_AlCl -> coeffs.dAlCl
            D_AlCl2 -> coeffs.dAlCl2
            D_AlCl3 -> coeffs.dAlCl3
            PG_H2 -> coeffs.pGH2
            PG_HCl -> coeffs.pGHCl
            PG_AlCl -> coeffs.pGAlCl
            PG_AlCl2 -> coeffs.pGAlCl2
            PG_AlCl3 -> coeffs.pGAlCl3
            else -> error(":(")
        }

    override val variableIdentifiers: List<PolinbCoefficient> = listOf(P_AlCl, P_H2, P_HCl, P_AlCl2, P_AlCl3)
}

