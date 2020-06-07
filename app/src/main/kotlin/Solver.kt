import kotlin.math.pow

interface Env<TOKEN> {
    /**
     * Returns the coefficient value by its identifier
     */
    val TOKEN.c: Double

    /**
     * List of available variable identifiers
     *
     * The order of variables, induced by this list, should coincide with parameter enumeration in the Jacobi matrix.
     */
    val variableIdentifiers: List<TOKEN>

    operator fun TOKEN.times(x: Double): Double = c * x
    operator fun Double.times(x: TOKEN): Double = this * x.c
    operator fun TOKEN.times(x: TOKEN): Double = c * x.c

    operator fun TOKEN.plus(x: Double): Double = c + x
    operator fun Double.plus(x: TOKEN): Double = this + x.c
    operator fun TOKEN.plus(x: TOKEN): Double = c + x.c

    operator fun TOKEN.minus(x: Double): Double = c - x
    operator fun Double.minus(x: TOKEN): Double = this - x.c
    operator fun TOKEN.minus(x: TOKEN): Double = c - x.c

    operator fun TOKEN.div(x: Double): Double = c / x
    operator fun Double.div(x: TOKEN): Double = this / x.c
    operator fun TOKEN.div(x: TOKEN): Double = c / x.c

    fun TOKEN.pow(x: Double): Double = c.pow(x)

    fun TOKEN.pow(x: Int): Double = c.pow(x)

    fun addIdentifiers(v: Vector): Env<TOKEN> = NestedEnv(v, this)
}

fun <T, U> intercept(f: Env<U>.(Vector) -> T): Env<U>.(Vector) -> T = { vec: Vector ->
    addIdentifiers(vec).f(vec)
}

fun <T, U> interceptedListOf(vararg fs: Env<U>.(Vector) -> T): List<Env<U>.(Vector) -> T> {
    return listOf(*fs).map(::intercept)
}

class NestedEnv<TOKEN>(private val additional: Vector, private val parentEnv: Env<TOKEN>) : Env<TOKEN> {
    override val TOKEN.c: Double
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

    override val variableIdentifiers: List<TOKEN>
        get() = parentEnv.variableIdentifiers

}

val zero: Env<*>.(Vector) -> Double = { 0.0 }


fun <U> Env<U>.generateSolution(
    system: Env<U>.(Vector) -> Vector,
    jacobi: Matrix<Env<U>.(Vector) -> Double>,
    startApproach: Vector
): Sequence<List<Double>> =
    generateSequence(startApproach) { x_k ->
        val fx = system(x_k)
        val matrix = jacobi.map { it.map { it(x_k) } }.inverted()
        val subtracted = matrix * fx
        x_k.zip(subtracted, Double::minus)
    }


