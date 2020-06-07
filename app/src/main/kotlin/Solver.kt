import kotlin.math.pow


interface Env<TOKEN, R> {

    val TOKEN.c: R
    val Double.c: R

    val variableIdentifiers: List<TOKEN>

    operator fun R.times(x: R): R

    operator fun R.plus(x: R): R

    operator fun R.minus(x: R): R

    operator fun R.div(x: R): R

    infix fun R.pow(x: Double): R

    infix fun R.pow(x: Int): R

    fun zero(): (List<*>) -> R = { 0.0.c }

    fun addIdentifiers(v: List<R>): Env<TOKEN, R> = NestedEnv(v, this)
}

fun Double.uniquePow(x: Double) = this.pow(x)

fun Double.uniquePow(x: Int) = this.pow(x)

@Suppress("EXTENSION_SHADOWED_BY_MEMBER")
interface EvaluationEnv<TOKEN> : Env<TOKEN, Double> {
    /**
     * Returns the coefficient value by its identifier
     */
    override val TOKEN.c: Double

    /**
     * List of available variable identifiers
     *
     * The order of variables, induced by this list, should coincide with parameter enumeration in the Jacobi matrix.
     */
    override val variableIdentifiers: List<TOKEN>
    override val Double.c: Double
        get() = this

    override fun Double.times(x: Double): Double = this.times(x)

    override fun Double.plus(x: Double): Double = this.plus(x)

    override fun Double.minus(x: Double): Double = this.minus(x)

    override fun Double.div(x: Double): Double = this.div(x)

    override fun Double.pow(x: Double): Double = this.uniquePow(x)

    override fun Double.pow(x: Int): Double = this.uniquePow(x)

}


fun <T, U, R> intercept(f: Env<U, R>.(List<R>) -> T): Env<U, R>.(List<R>) -> T = { vec: List<R> ->
    addIdentifiers(vec).f(vec)
}

fun <T, U, R> interceptedListOf(vararg fs: Env<U, R>.(List<R>) -> T): List<Env<U, R>.(List<R>) -> T> {
    return listOf(*fs).map(::intercept)
}

class NestedEnv<TOKEN, R>(private val additional: List<R>, private val parentEnv: Env<TOKEN, R>) :
    Env<TOKEN, R> by parentEnv {
    override val TOKEN.c: R
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

}


fun <U> generateSolution(
    env: EvaluationEnv<U>,
    system: Env<U, Any>.(List<Any>) -> List<Any>,
    jacobi: Matrix<Env<U, Any>.(List<Any>) -> Any>,
    startApproach: Vector
): Sequence<List<Double>> {

    return generateSequence(startApproach) { x_k: List<Double> ->
        @Suppress("UNCHECKED_CAST") val fx = (env as Env<U, Any>).system(x_k) as List<Double>
        @Suppress("UNCHECKED_CAST") val matrix =
            (jacobi.map { it.map { (env as Env<U, Any>).it(x_k) } } as Matrix<Double>).inverted()
        val subtracted: List<Double> = matrix * fx
        // here goes AST construction
        // val ASTEnv = ExpressionEnv(env)
        // val ast : Expression = ASTEnv.system(x_k - a * subtracted)
        // val a = ast.differentiate().findOptimum()
        x_k.zip(subtracted) { x, y -> x - /*a * */ y }
    }
}


