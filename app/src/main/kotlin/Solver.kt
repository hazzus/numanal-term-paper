interface Env {
    /**
     * Returns the coefficient value by its name
     */
    fun String.c(): Double

    /**
     * Returns variable value by provided name
     *
     * The order, generated by this function, should coincide with parameter enumeration in the Jacobian matrix.
     */
    operator fun List<Double>.get(key: String): Double

    operator fun String.times(x: Double): Double = c() * x
    operator fun Double.times(x: String): Double = this * x.c()

    operator fun String.minus(x: Double): Double = c() - x
    operator fun Double.minus(x: String): Double = this - x.c()
}

val zero: Env.(Vector) -> Double = { 0.0 }


fun <T> generateSolution(
    env: T,
    system: T.(Vector) -> Vector,
    jacobi: Matrix<T.(Vector) -> Double>,
    startApproach: Vector
): Sequence<List<Double>> =
    generateSequence(startApproach) { x_k ->
        val fx = env.system(x_k)
        val matrix = jacobi.map { it.map { env.it(x_k) } }.inverted()
        val subtracted = matrix * fx
        x_k.zip(subtracted, Double::minus)
    }


