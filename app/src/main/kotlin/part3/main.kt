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

fun solve(pGH2: Double, solutionWriter: OutputStreamWriter, plotDataDir: String) {
    val alCl3Writer = File(plotDataDir + "G_AlCl3_to_x_g").writer()
    val gaClWriter = File(plotDataDir + "G_GaCl_to_x_g").writer()
    val vWriter = File(plotDataDir + "V_g_AlGaN_to_x_g").writer()
    val xWriter = File(plotDataDir + "x_to_x_g").writer()

    var initialApproximation = listOf(0.3, 0.3, 0.3, 0.3, 0.3, 0.3)

    var pGAlCl3 = 0.0
    while (pGAlCl3 <= 30.0) {
        val coefficients = computeCoefficients(pGAlCl3, pGH2)

        val xG = pGAlCl3 / 30

        val env = Part3Env(coefficients)
        val eps = 1e-14
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

        // get x
        val x = solution.last()[0]

        // write results with description
//        solutionWriter.appendln("Solution for P^G_AlCl3 = $pGAlCl3, P^G_H2 = $pGH2:")
//        solutionWriter.appendln(solutionStr.substring(1, length - 1))
//        solutionWriter.appendln("G_AlCl3 = $gAlCl3, G_GaCl = $gGaCl")
//        solutionWriter.appendln("V^G_AlGaN = $vGAlGaN")
//        solutionWriter.appendln()
//
//        solutionWriter.flush()

        // write plot data

        alCl3Writer.appendln("${"%.10f".format(xG)} ${"%.15f".format(gAlCl3)}")
        gaClWriter.appendln("${"%.10f".format(xG)} ${"%.15f".format(gGaCl)}")
        vWriter.appendln("${"%.10f".format(xG)} ${"%.15f".format(vGAlGaN)}")
        xWriter.appendln("${"%.10f".format(xG)} ${"%.15f".format(x)}")

        alCl3Writer.flush()
        gaClWriter.flush()
        vWriter.flush()
        xWriter.flush()

        pGAlCl3 += 0.01
        initialApproximation = solution.last()
    }

    alCl3Writer.close()
    gaClWriter.close()
    vWriter.close()
    xWriter.close()
}

fun calculateResults() {
    val writer = File("app/part3Files/output/solution").writer()
    solve(0.0, writer, "app/part3Files/plots/firstCase/")
    solve(9847.0, writer, "app/part3Files/plots/secondCase/")
    writer.close()
}

fun main() {
    calculateResults()
}