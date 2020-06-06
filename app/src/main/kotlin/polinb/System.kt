package polinb

import Env
import Matrix
import Vector
import intercept
import interceptedListOf
import zero


const val K1 = "k1"
const val K2 = "k2"
const val K3 = "k3"
const val D_HCl = "d_hcl"
const val D_H2 = "d_h2"
const val D_AlCl = "d_alcl"
const val D_AlCl2 = "d_alcl2"
const val D_AlCl3 = "d_alcl3"
const val P_AlCl = "p_alcl"
const val P_H2 = "p_h2"
const val P_HCl = "p_hcl"
const val P_AlCl2 = "p_alcl2"
const val P_AlCl3 = "p_alcl3"
const val PG_AlCl = "pg_alcl"
const val PG_HCl = "pg_hcl"
const val PG_H2 = "pg_h2"
const val PG_AlCl2 = "pg_alcl2"
const val PG_AlCl3 = "pg_alcl3"

val polinbSystem: Env.(Vector) -> Vector = intercept {
    listOf(
        K1 * P_AlCl.pow(2.0) * P_H2 - P_HCl.pow(2.0),
        K2 * P_AlCl2 * P_H2 - P_HCl.pow(2.0),
        K3 * P_AlCl3.pow(2.0) * P_H2.pow(3.0) - P_HCl.pow(6),
        D_HCl * (PG_HCl - P_HCl) + 2.0 * D_H2 * (PG_H2 - P_H2),
        D_AlCl * (PG_AlCl - P_AlCl) + 2.0 * D_AlCl2 * (PG_AlCl2 - P_AlCl2) + 3.0 * D_AlCl3 * (PG_AlCl3 - P_AlCl3) + D_HCl * (PG_HCl - P_HCl)
    )
}

val polinbJacobi: Matrix<Env.(Vector) -> Double> = listOf(
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

class PolinbEnv(private val coeffs: Coefficients) : Env {
    override val String.c: Double
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

    override val variableIdentifiers: List<String> = listOf("p_alcl", "p_h2", "p_hcl", "p_alcl2", "p_alcl3")
}

