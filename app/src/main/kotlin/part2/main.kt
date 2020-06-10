package part2

import dot
import generateSolution
import kotlin.math.pow

class GResults(val gAlCl: Double,
               val gAlCl2: Double,
               val gAlCl3: Double) {
    fun print() = "G_AlCl = $gAlCl\n" +
            "G_AlCl2 = $gAlCl2\n" +
            "G_AlCl3 = $gAlCl3"
}

class PResults(val pAlCl: Double,
               val pAlCl2: Double,
               val pAlCl3: Double) {
    fun print() = "P_AlCl = $pAlCl\n" +
            "P_AlCl2 = $pAlCl2\n" +
            "P_AlCl3 = $pAlCl3"
}

fun computeG(T: Double, coefficients: Coefficients, ps: PResults) : GResults {
    val r = 8314.0
    val d = 0.01
    val denominator = r * T * d
    val gAlCl = coefficients.dGaCl * (coefficients.pGGaCl - ps.pAlCl) / denominator
    val gAlCl2 = coefficients.dGaCl2 * (coefficients.pGGaCl2 - ps.pAlCl2) / denominator
    val gAlCl3 = coefficients.dGaCl3 * (coefficients.pGGaCl3 - ps.pAlCl3) / denominator
    return GResults(gAlCl, gAlCl2, gAlCl3)
}

fun computeV(T: Double): Double {
    val res = computeCoefficients(T)
//    println(res.print())
    // solving a system of equations
    val env = Part2Env(res)
    val startApproach = listOf(0.1, 0.1, 0.1, 0.1, 0.1)
    val eps = 1e-20
    val vectorSequence = env.generateSolution(
        part2System,
        part2Jacobi,
        startApproach
    )
    val solution =
        vectorSequence.takeWhile { vec -> env.part2System(vec).run { dot(this, this) } > eps }.toList()
    val pRes = PResults(solution.last()[0], solution.last()[3], solution.last()[4])
//    println(pRes.print())
//    println()
//    solution.last().zip(env.variableIdentifiers).forEach { (value, name) ->
//        println("$name: $value")
//    }
    val gRes = computeG(T, res, pRes)
    println(gRes.print())
    println()
    val mu = 26.9815 // kg/kmol
    val po = 2690.0 // kg/m^3
    return (gRes.gAlCl + gRes.gAlCl2 + gRes.gAlCl3) * mu * 10.0.pow(9) / po
}

fun main() {
    val t = 623.15
    println("V = ${computeV(t)}")
}
