package polinb

import Env
import EvaluationEnv
import Matrix
import intercept
import interceptedListOf
import polinb.PolinbCoefficient.*

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

fun <R> polinbSystem(): Env<PolinbCoefficient, R>.(List<R>) -> List<R> = intercept {
    listOf(
        K1.c * (P_AlCl.c.pow(2)) * P_H2.c - (P_HCl.c pow 2),
        K2.c * P_AlCl2.c * P_H2.c - (P_HCl.c pow 2),
        K3.c * (P_AlCl3.c pow 2) * (P_H2.c pow 3) - (P_HCl.c pow 6),
        D_HCl.c * (PG_HCl.c - P_HCl.c) + 2.0.c * D_H2.c * (PG_H2.c - P_H2.c),
        D_AlCl.c * (PG_AlCl.c - P_AlCl.c) + 2.0.c * D_AlCl2.c * (PG_AlCl2.c - P_AlCl2.c) + 3.0.c * D_AlCl3.c * (PG_AlCl3.c - P_AlCl3.c) + D_HCl.c * (PG_HCl.c - P_HCl.c)
    )
}

fun <R> polinbJacobi(): Matrix<Env<PolinbCoefficient, R>.(List<R>) -> R> = listOf(
    interceptedListOf(
        { 2.0.c * K1.c * P_AlCl.c * P_H2.c }, // alcl
        { K1.c * P_AlCl.c.pow(2.0) }, // h2
        { (-2.0).c * P_HCl.c }, // hcl
        { 0.0.c }, // alcl2
        { 0.0.c } // alcl3
    ),
    interceptedListOf(
        { 0.0.c },
        { K2.c * P_AlCl2.c },
        { (-2.0).c * P_HCl.c },
        { K2.c * P_H2.c },
        { 0.0.c }
    ),
    interceptedListOf(
        { 0.0.c },
        { 3.0.c * K3.c * P_AlCl3.c.pow(2.0) * P_H2.c.pow(2.0) },
        { (-6.0).c * P_HCl.c.pow(5.0) },
        { 0.0.c },
        { K3.c * 2.0.c * P_AlCl3.c * P_H2.c.pow(3.0) }),
    interceptedListOf(
        { 0.0.c },
        { (-2.0).c * D_H2.c },
        { D_HCl.c * (-1.0).c },
        { 0.0.c },
        { 0.0.c }
    ),
    interceptedListOf(
        { D_AlCl.c * (-1.0).c },
        { 0.0.c },
        { (-1.0).c * D_HCl.c },
        { (-2.0).c * D_AlCl2.c },
        { (-3.0).c * D_AlCl3.c }
    )
)

class PolinbEnv(private val coeffs: Coefficients) : EvaluationEnv<PolinbCoefficient> {
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

