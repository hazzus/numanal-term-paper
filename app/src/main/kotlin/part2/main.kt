package part2

import dot
import generateSolution
import kotlin.math.pow

class GResults(val gGaCl: Double,
               val gGaCl2: Double,
               val gGaCl3: Double) {
    fun print() = "G_GaCl = $gGaCl\n" +
            "G_GaCl2 = $gGaCl2\n" +
            "G_GaCl3 = $gGaCl3"
}

class PResults(val pGaCl: Double,
               val pGaCl2: Double,
               val pGaCl3: Double) {
    fun print() = "P_GaCl = $pGaCl\n" +
            "P_GaCl2 = $pGaCl2\n" +
            "P_GaCl3 = $pGaCl3"
}

fun computeG(T: Double, coefficients: Coefficients, ps: PResults) : GResults {
    val r = 8314.0
    val d = 0.01
    val denominator = r * T * d
    val gGaCl = coefficients.dGaCl * (coefficients.pGGaCl - ps.pGaCl) / denominator
    val gGaCl2 = coefficients.dGaCl2 * (coefficients.pGGaCl2 - ps.pGaCl2) / denominator
    val gGaCl3 = coefficients.dGaCl3 * (coefficients.pGGaCl3 - ps.pGaCl3) / denominator
    return GResults(gGaCl, gGaCl2, gGaCl3)
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
    return (gRes.gGaCl + gRes.gGaCl2 + gRes.gGaCl3) * mu * 10.0.pow(9) / po
}

fun main() {
    val t = 923.15
    println("V = ${computeV(t)}")
}
