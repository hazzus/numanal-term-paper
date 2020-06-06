val polinbSystem: Env.(Vector) -> Vector = intercept {
    listOf(
        "k1" * "p_alcl".pow(2.0) * "p_h2" - "p_hcl".pow(2.0),
        "k2" * "p_alcl2" * "p_h2" - "p_hcl".pow(2.0),
        "k3" * "p_alcl3".pow(2.0) * "p_h2".pow(3.0) - "p_hcl".pow(6),
        "d_hcl" * ("pg_hcl" - "p_hcl") + 2.0 * "d_h2" * ("pg_h2" - "p_h2"),
        "d_alcl" * ("pg_alcl" - "p_alcl") + 2.0 * "d_alcl2" * ("pg_alcl2" - "p_alcl2") + 3.0 * "d_alcl3" * ("pg_alcl3" - "p_alcl3") + "d_hcl" * ("pg_hcl" - "p_hcl")
    )
}

val polinbJacobi: Matrix<Env.(Vector) -> Double> = listOf(
    interceptedListOf(
        { 2.0 * "k1" * "p_alcl" * "p_h2" }, // alcl
        { "k1" * "p_alcl".pow(2.0) }, // h2
        { -2.0 * "p_hcl" }, // hcl
        zero, // alcl2
        zero // alcl3
    ),
    interceptedListOf(
        zero,
        { "k2" * "p_alcl2" },
        { -2.0 * "p_hcl" },
        { "k2" * "p_h2" },
        zero
    ),
    interceptedListOf(
        zero,
        { 3.0 * "k3" * "p_alcl3".pow(2.0) * "p_h2".pow(2.0) },
        { -6.0 * "p_hcl".pow(5.0) },
        zero,
        { "k3" * 2.0 * "p_alcl3" * "p_h2".pow(3.0) }),
    interceptedListOf(
        zero,
        { -2.0 * "d_h2" },
        { "d_hcl" * (-1.0) },
        zero,
        zero
    ),
    interceptedListOf(
        { "d_alcl" * (-1.0) },
        zero,
        { (-1.0) * "d_hcl" },
        { -2.0 * "d_alcl2" },
        { -3.0 * "d_alcl3" }
    )
)

object PolinbEnv : Env {
    override val String.c: Double
        get() {
            return 0.001
        }

    override val variableIdentifiers: List<String> = listOf("p_alcl", "p_h2", "p_hcl", "p_alcl2", "p_alcl3")
}


fun main() {
    val startApproach = listOf(-0.1, -0.1, -0.1, -0.1, -0.1)
    val eps = 1e-6
    val vectorSequence = generateSolution(PolinbEnv, polinbSystem, polinbJacobi, startApproach)
    val solution =
        vectorSequence.takeWhile { vec -> PolinbEnv.polinbSystem(vec).run { dot(this, this) } > eps }.toList()
    println(solution)
}