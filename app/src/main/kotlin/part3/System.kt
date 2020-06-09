package part3

import Env
import Matrix
import Vector
import intercept
import interceptedListOf
import part3.Part3Coefficient.*
import part3.coefficients.Coefficients
import zero

@Suppress("EnumEntryName")
enum class Part3Coefficient {
    K9,
    K10,
    D_HCl,
    D_H2,
    D_NH3,
    D_AlCl3,
    D_GaCl,
    PE_AlCl3,
    PE_NH3,
    PE_HCl,
    PE_H2,
    PE_GaCl,
    PG_AlCl3,
    PG_NH3,
    PG_HCl,
    PG_H2,
    PG_GaCl,
    x
}

val part3System: Env<Part3Coefficient>.(Vector) -> Vector = intercept {
    listOf(
            K9 * x * PE_HCl.pow(3.0) - PE_AlCl3 * PE_NH3,
            K10 * (1.0 - x) * PE_HCl * PE_H2 - PE_GaCl * PE_NH3,
            D_HCl * (PG_HCl - PE_HCl) + 2.0 * D_H2 * (PG_H2 - PE_H2) + 3.0 * D_NH3 * (PG_NH3 - PE_NH3),
            3.0 * D_AlCl3 * (PG_AlCl3 - PE_AlCl3) + D_GaCl * (PG_GaCl - PE_GaCl) + D_HCl * (PG_HCl - PE_HCl),
            D_AlCl3 * (PG_AlCl3 - PE_AlCl3) + D_GaCl * (PG_GaCl - PE_GaCl) - D_NH3 * (PG_NH3 - PE_NH3),
            D_AlCl3 * (PG_AlCl3 - PE_AlCl3) * (1.0 - x) - D_GaCl * (PG_GaCl - PE_GaCl) * x
    )
}

val part3Jacobi: Matrix<Env<Part3Coefficient>.(Vector) -> Double> = listOf(
        interceptedListOf(
                { K9 * PE_HCl.pow(3.0) },
                { -1.0 * PE_NH3 },
                zero,
                { -1.0 * PE_AlCl3 },
                { 3.0 * K9 * x * (PE_HCl).pow(2.0) },
                zero
        ),
        interceptedListOf(
                { -1.0 * K10 * PE_HCl * PE_H2 },
                zero,
                { -1.0 * PE_NH3 },
                { -1.0 * PE_GaCl },
                { K10 * (1.0 - x) * PE_H2 },
                { K10 * (1.0 - x) * PE_HCl }
        ),
        interceptedListOf(
                zero,
                zero,
                zero,
                { -3.0 * D_NH3 },
                { -1.0 * D_HCl },
                { -2.0 * D_H2 }
        ),
        interceptedListOf(
                zero,
                { -3.0 * D_AlCl3 },
                { -1.0 * D_GaCl },
                zero,
                { -1.0 * D_HCl },
                zero
        ),
        interceptedListOf(
                zero,
                { -1.0 * D_AlCl3 },
                { -1.0 * D_GaCl },
                { 1.0 * D_NH3 },
                zero,
                zero
        ),
        interceptedListOf(
                { -1.0 * D_AlCl3 * (PG_AlCl3 - PE_AlCl3) - D_GaCl * (PG_GaCl - PE_GaCl) },
                { (x - 1.0) * D_AlCl3 },
                { D_GaCl * x },
                zero,
                zero,
                zero
        )
)

class Part3Env(private val coeffs: Coefficients) : Env<Part3Coefficient> {
    override val Part3Coefficient.c: Double
        get() = when (this) {
            K9 -> coeffs.K9
            K10 -> coeffs.K10
            D_HCl -> coeffs.dHCl
            D_H2 -> coeffs.dH2
            D_NH3 -> coeffs.dNH3
            D_AlCl3 -> coeffs.dAlCl3
            D_GaCl -> coeffs.dGaCl
            PG_AlCl3 -> coeffs.pGAlCl3
            PG_H2 -> coeffs.pGH2
            PG_GaCl -> coeffs.pGGaCl
            PG_HCl -> coeffs.pGHCl
            PG_NH3 -> coeffs.pGNH3
            else -> error(":(")
        }

    override val variableIdentifiers: List<Part3Coefficient> = listOf(x, PE_AlCl3, PE_GaCl, PE_NH3, PE_HCl, PE_H2)
}

