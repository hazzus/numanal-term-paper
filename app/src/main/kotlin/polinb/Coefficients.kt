package polinb

import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sqrt


class Coefficients(
    val dAlCl: Double,
    val dAlCl2: Double,
    val dAlCl3: Double,
    val dHCl: Double,
    val dH2: Double,
    val K1: Double,
    val K2: Double,
    val K3: Double,
    val pGAlCl: Double = 0.0,
    val pGAlCl2: Double = 0.0,
    val pGAlCl3: Double = 0.0,
    val pGH2: Double = 0.0,
    val pGHCl: Double = 10_000.0,
    val pGN2: Double = 90_000.0
) {
    fun print() = "D_AlCl = $dAlCl\n" +
            "D_AlCl2 = $dAlCl2\n" +
            "D_AlCl3 = $dAlCl3\n" +
            "D_HCl = $dHCl\n" +
            "D_H2 = $dH2\n" +
            "K1 = $K1\n" +
            "K2 = $K2\n" +
            "K3 = $K3\n"
}

fun computeCoefficients(T: Double): Coefficients {
    val elAlCl = AlCl()
    val elAlCl2 = AlCl2()
    val elAlCl3 = AlCl3()
    val elHCl = HCl()
    val elH2 = H2()
    val elAl = Al()
    val gAlCl = elAlCl.computeG0(T)
    val gAlCl2 = elAlCl2.computeG0(T)
    val gAlCl3 = elAlCl3.computeG0(T)
    val gH2 = elH2.computeG0(T)
    val gAl = elAl.computeG0(T)
    val gHCl = elHCl.computeG0(T)
    val r = 8.314
    val p = 100_000.0
    val k1 = exp((2 * gAlCl + gH2 - 2 * gAl - 2 * gHCl) / (r * T)) / p
    val k2 = exp((gAlCl2 + gH2 - gAl - 2 * gHCl) / (r * T))
    val k3 = exp((2 * gAlCl3 + 3 * gH2 - 2 * gAl - 6 * gHCl) / (r * T)) * p
    return Coefficients(
        elAlCl.computeD(T),
        elAlCl2.computeD(T),
        elAlCl3.computeD(T),
        elHCl.computeD(T),
        elH2.computeD(T),
        k1, k2, k3
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

class AlCl : Element() {
    override val sigma = 3.58
    override val epsilon = 932.0
    override val mu = 62.4345
    override val name = "AlCl"
    override val h = -51032.0
    override val f1 = 318.9948
    override val f2 = 36.94626
    override val f3 = -0.001226431
    override val f4 = 1.1881743
    override val f5 = 5.638541
    override val f6 = -5.066135
    override val f7 = 5.219347
}

class AlCl2 : Element() {
    override val sigma = 5.3
    override val epsilon = 825.0
    override val mu = 97.8875
    override val name = "AlCl2"
    override val h = -259000.0
    override val f1 = 427.2137
    override val f2 = 56.56409
    override val f3 = -0.002961273
    override val f4 = 1.893842
    override val f5 = 12.40072
    override val f6 = -22.65441
    override val f7 = 21.29898
}

class AlCl3 : Element() {
    override val sigma = 5.13
    override val epsilon = 472.0
    override val mu = 133.3405
    override val name = "AlCl3"
    override val h = -584100.0
    override val f1 = 511.8114
    override val f2 = 81.15042
    override val f3 = -0.004834879
    override val f4 = 2.752097
    override val f5 = 13.40078
    override val f6 = -21.28001
    override val f7 = 16.92868
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

class Al : Element() {
    override val sigma = 0.0
    override val epsilon = 0.0
    override val mu = 0.0
    override val name = "Al"
    override val h = 0.0
    override val f1 = 172.8289
    override val f2 = 50.51806
    override val f3 = -0.00411847
    override val f4 = 1.476107
    override val f5 = -458.1279
    override val f6 = 2105.75
    override val f7 = -4168.337
}

