package part2

import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sqrt


class Coefficients(
        val dHCl: Double,
        val dH2: Double,
        val dGaCl: Double,
        val dGaCl2: Double,
        val dGaCl3: Double,
        val K4: Double,
        val K5: Double,
        val K6: Double,
        // P^g_i
        val pGGaCl: Double = 0.0,
        val pGGaCl2: Double = 0.0,
        val pGGaCl3: Double = 0.0,
        val pGH2: Double = 0.0,
        val pGHCl: Double = 10_000.0,
        val pGN2: Double = 90_000.0
) {
    fun print() = "D_HCl = $dHCl\n" +
            "D_H2 = $dH2\n" +
            "D_GaCl = $dGaCl\n" +
            "D_GaGl2 = $dGaCl2\n" +
            "D_GaCl3 = $dGaCl3\n" +
            "K4 = $K4\n" +
            "K5 = $K5\n" +
            "K6 = $K6\n"
}

fun computeCoefficients(T: Double): Coefficients {
    val elGaCl = GaCl()
    val elGaCl2 = GaCl2()
    val elGaCl3 = GaCl3()
    val elHCl = part2.HCl()
    val elH2 = part2.H2()
    val elGa = Ga()
    val gGaCl = elGaCl.computeG0(T)
    val gGaCl2 = elGaCl2.computeG0(T)
    val gGaCl3 = elGaCl3.computeG0(T)
    val gH2 = elH2.computeG0(T)
    val gAl = elGa.computeG0(T)
    val gHCl = elHCl.computeG0(T)
    val r = 8.314
    val p = 100_000.0
    // val k1 = exp((2 * gAlCl + gH2 - 2 * gAl - 2 * gHCl) / (r * T)) / p
    // val k2 = exp((gAlCl2 + gH2 - gAl - 2 * gHCl) / (r * T))
    // val k3 = exp((2 * gAlCl3 + 3 * gH2 - 2 * gAl - 6 * gHCl) / (r * T)) * p
    val k4 = 1.0
    val k5 = 1.0
    val k6 = 1.0
    return part2.Coefficients(
            elGaCl.computeD(T),
            elGaCl2.computeD(T),
            elGaCl3.computeD(T),
            elHCl.computeD(T),
            elH2.computeD(T),
            k4, k5, k6
    )
}

abstract class Element {
    abstract val sigma: Double
    abstract val epsilon: Double
    abstract val mu: Double
    abstract val name: String
    abstract val h: Double
    abstract val f1: Double
    abstract val f2: Double
    abstract val f3: Double
    abstract val f4: Double
    abstract val f5: Double
    abstract val f6: Double
    abstract val f7: Double

    fun computeD(T: Double): Double {
        val sigmaN2 = 3.798
        val p = 100000
        val sigmaElN2 = (sigma + sigmaN2) / 2
        val epsilonN2 = 71.4
        val epsilonElN2 = sqrt(epsilonN2 * epsilon)
        val omega = 1.074 * (T / epsilonElN2).pow(-0.1604)
        val muN2 = 28.0135
        val muElN2 = 2 * muN2 * mu / (muN2 + mu)
        return 2.628 * 10.0.pow(-2.0) * T.pow(1.5) / (p * sigmaElN2 * omega * muElN2.pow(0.5))
    }

    fun printD(T: Double): String {
        val d = computeD(T)
        return "D_$name = $d"
    }

    fun computeF(T: Double): Double {
        val x = T / (10.0.pow(4.0))
        return f1 + f2 * ln(x) + f3 * x.pow(-2.0) + f4 * (1.0 / x) + f5 * x + f6 * x * x + f7 * x * x * x
    }

    fun computeG0(T: Double): Double {
        return h - computeF(T) * T
    }

    fun printG0(T: Double): String {
        val g = computeG0(T)
        return "G0_$name = $g"
    }
}

class HCl : Element() {
    override val sigma = 2.737
    override val epsilon = 167.1
    override val mu = 36.461
    override val name = "HCl"
    override val h = -92310.0
    override val f1 = 243.9878
    override val f2 = 23.15984
    override val f3 = 0.001819985
    override val f4 = 0.6147384
    override val f5 = 51.16604
    override val f6 = -36.89502
    override val f7 = 9.174252
}

class H2 : Element() {
    override val sigma = 2.93
    override val epsilon = 34.1
    override val mu = 2.016
    override val name = "H2"
    override val h = 0.0
    override val f1 = 205.5368
    override val f2 = 29.50487
    override val f3 = 0.000168424
    override val f4 = 0.86065612
    override val f5 = -14.95312
    override val f6 = 78.18955
    override val f7 = -82.78981
}

class GaCl : Element() {
    override val sigma = 3.696
    override val epsilon = 348.2
    override val mu = 105.173
    override val name = "GaCl"
    override val h = -70553.0
    override val f1 = 332.2718
    override val f2 = 37.11052
    override val f3 = -0.000746187
    override val f4 = 1.1606512
    override val f5 = 4.891346
    override val f6 = -4.467591
    override val f7 = 5.506236
}

class GaCl2 : Element() {
    override val sigma = 4.293
    override val epsilon = 465.0
    override val mu = 140.626
    override val name = "GaCl2"
    override val h = -241238.0
    override val f1 = 443.2976
    override val f2 = 57.745845
    override val f3 = -0.002265112
    override val f4 = 1.8755545
    override val f5 = 3.66186
    override val f6 = -9.356338
    override val f7 = 15.88245
}

class GaCl3 : Element() {
    override val sigma = 5.034
    override val epsilon = 548.24
    override val mu = 176.080
    override val name = "GaCl3"
    override val h = -431573.0
    override val f1 = 526.8113
    override val f2 = 82.03355
    override val f3 = -0.003486473
    override val f4 = 2.6855923
    override val f5 = 8.278878
    override val f6 = -14.5678
    override val f7 = 12.8899
}

class Ga : Element() {
    override val sigma = 69.723
    override val epsilon = 0.0
    override val mu = 0.0
    override val name = "Ga"
    override val h = 0.0
    override val f1 = 125.9597
    override val f2 = 26.03107
    override val f3 = 0.001178297
    override val f4 = 0.13976
    override val f5 = -0.5698425
    override val f6 = 0.04723008
    override val f7 = 7.212525
}