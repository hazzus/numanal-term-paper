package part3.coefficients

import java.io.File
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sqrt

const val P = 100_000
const val R = 8.314
const val T = 1373.15
const val DELIMITER = ';'
const val WHITESPACE = ' '

const val F_CONSTANTS = "app/part3Files/resources/f_constants"
const val H_CONSTANTS = "app/part3Files/resources/h_constants"
const val SIGMA_CONSTANTS = "app/part3Files/resources/sigma_constants"
const val EPSILON_CONSTANTS = "app/part3Files/resources/epsilon_constants"
const val MU_CONSTANTS = "app/part3Files/resources/mu_constants"

val variablesForDiffusion = listOf("AlCl3", "GaCl", "NH3", "H2", "HCl")
val variablesForPhiAndG = listOf("AlCl3", "GaCl", "NH3", "H2", "HCl", "AlN", "GaN")

val h = HashMap<String, Double>()
val sigma = HashMap<String, Double>()
val epsilon = HashMap<String, Double>()
val mu = HashMap<String, Double>()
val f = ArrayList<MutableMap<String, Double>>()

var diffusionCoefficients = HashMap<String, Double>()
var phi = HashMap<String, Double>()
var G = HashMap<String, Double>()

fun readFConstants() {
    File(F_CONSTANTS).readLines().forEach {
        val map = HashMap<String, Double>()
        for (pair in it.split(DELIMITER)) {
            val whitespaceIndex = pair.indexOf(WHITESPACE)
            val key = pair.substring(0, whitespaceIndex)
            val value = pair.substring(whitespaceIndex + 1).toDouble()
            map[key] = value
        }
        f.add(map)
    }
}

fun readConstants(fileName: String, map: MutableMap<String, Double>) {
    File(fileName).readLines()[0].split(DELIMITER).forEach {
        val whitespaceIndex = it.indexOf(WHITESPACE)
        val key = it.substring(0, whitespaceIndex)
        val value = it.substring(whitespaceIndex + 1).toDouble()
        map[key] = value
    }
}

fun initializeConstants() {
    readFConstants()
    readConstants(H_CONSTANTS, h)
    readConstants(SIGMA_CONSTANTS, sigma)
    readConstants(EPSILON_CONSTANTS, epsilon)
    readConstants(MU_CONSTANTS, mu)
}

fun calculatePhi(t: Double, i: String): Double {
    val x = t / 10_000
    return f[0][i]!! + f[1][i]!! * ln(x) + f[2][i]!! * x.pow(-2) +
            f[3][i]!! * x.pow(-1) + f[4][i]!! * x +
            f[5][i]!! * x.pow(2) + f[6][i]!! * x.pow(3)
}

fun calculateG(t: Double, i: String): Double = h[i]!! - phi[i]!! * t

fun calculateDiffusionCoefficient(t: Double, i: String): Double =
        0.02628 * t.pow(1.5) / (P * ((sigma[i]!! + sigma["N2"]!!) / 2) *
                1.074 * (t / sqrt(epsilon[i]!! * epsilon["N2"]!!)).pow(-0.1604)
                * sqrt(2 * mu[i]!! * mu["N2"]!! / (mu[i]!! + mu["N2"]!!)))

fun calculateDiffusionCoefficients(t: Double) =
        variablesForDiffusion.forEach { diffusionCoefficients[it] = calculateDiffusionCoefficient(t, it) }

fun calculatePhi(t: Double) = variablesForPhiAndG.forEach { phi[it] = calculatePhi(t, it) }

fun calculateG(t: Double) = variablesForPhiAndG.forEach { G[it] = calculateG(t, it) }

fun calculateConstants(t: Double) {
    calculateDiffusionCoefficients(t)
    calculatePhi(t)
    calculateG(t)
}

fun deltaG9(): Double = G["AlCl3"]!! + G["NH3"]!! - G["AlN"]!! - 3 * G["HCl"]!!

fun deltaG10(): Double = G["GaCl"]!! + G["NH3"]!! - G["GaN"]!! - G["HCl"]!! - G["H2"]!!

fun computeCoefficients(pGAlCl3: Double, pGH2: Double): Coefficients {
    initializeConstants()
    calculateConstants(T)
    val k9 = exp(-deltaG9() / (R * T)) / P
    val k10 = exp(-deltaG10() / (R * T))

    return Coefficients(
            diffusionCoefficients["HCl"]!!,
            diffusionCoefficients["H2"]!!,
            diffusionCoefficients["NH3"]!!,
            diffusionCoefficients["AlCl3"]!!,
            diffusionCoefficients["GaCl"]!!,
            k9,
            k10,
            pGAlCl3,
            pGH2,
            98470.0 - pGH2,
            30.0 - pGAlCl3
    )
}

class Coefficients(
        val dHCl: Double,
        val dH2: Double,
        val dNH3: Double,
        val dAlCl3: Double,
        val dGaCl: Double,
        val K9: Double,
        val K10: Double,
        val pGAlCl3: Double,
        val pGH2: Double,
        val pGN2: Double,
        val pGGaCl: Double,
        val pGHCl: Double = 0.0,
        val pGNH3: Double = 1500.0
) {
    fun print() = """
            dHCl = $dHCl
            dH2 = $dH2
            dNH3 = $dNH3
            dAlCl3 = $dAlCl3
            dGaCl = $dGaCl
            K9 = $K9
            K10 = $K10
            pGHCl = $pGHCl
            pGH2 = $pGH2
            pGN2 = $pGN2
            pGNH3 = $pGNH3
            pGAlCl3 = $pGAlCl3
            pGGaCl = $pGGaCl
            """.trimIndent()
}

