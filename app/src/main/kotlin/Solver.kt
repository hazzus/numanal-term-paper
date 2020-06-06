import kotlin.math.pow

interface Env {
    /**
     * Returns the coefficient value by its name
     */
    val String.c: Double

    /**
     * List of available variable names
     *
     * The order of variables, induced by this list, should coincide with parameter enumeration in the Jacobi matrix.
     */
    val variableIdentifiers: List<String>

    operator fun List<Double>.get(key: String): Double = get(variableIdentifiers.indexOf(key))

    operator fun String.times(x: Double): Double = c * x
    operator fun Double.times(x: String): Double = this * x.c
    operator fun String.times(x: String): Double = c * x.c

    operator fun String.plus(x: Double): Double = c + x
    operator fun Double.plus(x: String): Double = this + x.c
    operator fun String.plus(x: String): Double = c + x.c

    operator fun String.minus(x: Double): Double = c - x
    operator fun Double.minus(x: String): Double = this - x.c
    operator fun String.minus(x: String): Double = c - x.c

    fun String.pow(x: Double): Double = c.pow(x)

    fun String.pow(x: Int): Double = c.pow(x)

    fun addIdentifiers(v: Vector): Env = NestedEnv(v, this)
}

fun <T> intercept(f: Env.(Vector) -> T): Env.(Vector) -> T = { vec: Vector ->
    addIdentifiers(vec).f(vec)
}

fun <T> interceptedListOf(vararg fs: Env.(Vector) -> T): List<Env.(Vector) -> T> {
    return listOf(*fs).map(::intercept)
}

class NestedEnv(private val additional: Vector, private val parentEnv: Env) : Env {
    override val String.c: Double
        get() {
            val key = this
            return if (key in parentEnv.variableIdentifiers) {
                additional[parentEnv.variableIdentifiers.indexOf(key)]
            } else {
                with(parentEnv) {
                    key.c
                }
            }
        }

    override val variableIdentifiers: List<String>
        get() = parentEnv.variableIdentifiers

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


