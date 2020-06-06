package polinb

import dot
import generateSolution

fun main() {
    val t = 623.15
    val res = computeCoefficients(t)
    val env = PolinbEnv(res)
    val startApproach = listOf(-0.1, -0.1, -0.1, -0.1, -0.1)
    val eps = 1e-6
    val vectorSequence = generateSolution(
        env,
        polinbSystem,
        polinbJacobi,
        startApproach
    )
    val solution =
        vectorSequence.takeWhile { vec -> env.polinbSystem(vec).run { dot(this, this) } > eps }.toList()
    solution.last().zip(env.variableIdentifiers).forEach { (value, name) ->
        println("$name: $value")
    }

    println(res.print())
}
