import kotlin.math.pow

val polinbSystem: Env.(Vector) -> Vector = { args ->
    listOf(
        "k1" * args["p_alcl"].pow(2.0) * args["p_h2"] - args["p_hcl"].pow(2.0),
        "k2" * args["p_alcl2"] * args["p_h2"] - args["p_hcl"].pow(2.0),
        "k3" * args["p_alcl3"].pow(2.0) * args["p_h2"].pow(3.0) - args["p_hcl"].pow(6),
        "d_hcl" * ("pg_hcl" - args["p_hcl"]) + 2.0 * "d_h2" * ("pg_h2" - args["p_h2"]),
        "d_alcl" * ("pg_alcl" - args["p_alcl"]) + 2.0 * "d_alcl2" * ("pg_alcl2" - args["p_alcl2"]) + 3.0 * "d_alcl3" * ("pg_alcl3" - args["p_alcl3"]) + "d_hcl" * ("pg_hcl" - args["p_hcl"])
    )
}

val polinbJacobi: Matrix<Env.(Vector) -> Double> = listOf(
    listOf(
        { 2.0 * "k1" * it["p_alcl"] * it["p_h2"] }, // alcl
        { "k1" * it["p_alcl"].pow(2.0) }, // h2
        { -2.0 * it["p_hcl"] }, // hcl
        zero, // alcl2
        zero // alcl3
    ),
    listOf(
        zero,
        { "k2" * it["p_alcl2"] },
        { -2.0 * it["p_hcl"] },
        { "k2" * it["p_h2"] },
        zero
    ),
    listOf(
        zero,
        { 3.0 * "k3" * it["p_alcl3"].pow(2.0) * it["p_h2"].pow(2.0) },
        { -6.0 * it["p_hcl"].pow(5.0) },
        zero,
        { "k3" * 2.0 * it["p_alcl3"] * it["p_h2"].pow(3.0) }),
    listOf(
        zero,
        { -2.0 * "d_h2" },
        { "d_hcl" * (-1.0) },
        zero,
        zero
    ),
    listOf(
        { "d_alcl" * (-1.0) },
        zero,
        { (-1.0) * "d_hcl" },
        { -2.0 * "d_alcl2" },
        { -3.0 * "d_alcl3" }
    )
)

object PolinbEnv : Env {
    override fun String.c(): Double {
        return 0.001
    }

    private val variableOrder = listOf("p_alcl", "p_h2", "p_hcl", "p_alcl2", "p_alcl3")

    override fun List<Double>.get(key: String): Double = get(variableOrder.indexOf(key))
}


fun main() {
    val startApproach = listOf(-0.1, -0.1, -0.1, -0.1, -0.1)
    val eps = 1e-6
    val vectorSequence = generateSolution(PolinbEnv, polinbSystem, polinbJacobi, startApproach)
    val solution =
        vectorSequence.takeWhile { vec -> PolinbEnv.polinbSystem(vec).run { dot(this, this) } > eps }.toList()
    println(solution)
}