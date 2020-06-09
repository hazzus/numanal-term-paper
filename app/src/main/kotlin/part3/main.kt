package part3

import dot
import generateSolution
import part3.coefficients.T
import part3.coefficients.computeCoefficients
import part3.coefficients.diffusionCoefficients
import java.io.File
import java.io.OutputStreamWriter
import kotlin.math.pow

const val R = 8314.0
const val delta = 0.01
const val muAlN = 40.988
const val roAlN = 3200.0
const val muGaN = 83.730
const val roGaN = 6150.0

fun solve(pGH2: Double, writer: OutputStreamWriter) {
    var pGAlCl3 = 0.0
    while (pGAlCl3 <= 30.0) {
        println(pGAlCl3)
        val coefficients = computeCoefficients(pGAlCl3, pGH2)
        val env = Part3Env(coefficients)
        val initialApproximation = listOf(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
        val eps = 1e-15
        val vectorSequence = env.generateSolution(
                part3System,
                part3Jacobi,
                initialApproximation
        )

        val solution = vectorSequence.takeWhile { vec -> env.part3System(vec).run { dot(this, this) } > eps }.toList()
        val solutionStr = solution.last().toString()
        val length = solutionStr.length

        // calculate G_AlCl3 and G_GaCl
        val gAlCl3 = diffusionCoefficients["AlCl3"]!! * (coefficients.pGAlCl3 - solution.last()[1]) / (R * T * delta)
        val gGaCl = diffusionCoefficients["GaCl"]!! * (coefficients.pGGaCl - solution.last()[2]) / (R * T * delta)

        // calculate V^G_AlGaN
        val vGAlGaN = (gAlCl3 * (muAlN / roAlN) + gGaCl * (muGaN / roGaN)) * 10.0.pow(9)

        writer.appendln("Solution for P^G_AlCl3 = $pGAlCl3, P^G_H2 = $pGH2:")
        writer.appendln(solutionStr.substring(1, length - 1))
        writer.appendln("G_AlCl3 = $gAlCl3, G_GaCl = $gGaCl")
        writer.appendln("V^G_AlGaN = $vGAlGaN")
        writer.appendln()

        writer.flush()

        pGAlCl3 += 0.001
    }
}

fun calculateResults() {
    val writer = File("app/part3Files/output/solution").writer()
    solve(0.0, writer)
    solve(9847.0, writer)
    writer.close()
}

fun main() {
    calculateResults()
}