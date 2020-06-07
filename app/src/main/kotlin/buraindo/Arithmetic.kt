package buraindo

import kotlin.math.pow

fun varOf(s: String): Var = Var(s)
fun const(n: Double): Const = Const(n)
fun expression(e: Any): Expression = when (e) {
    is Double -> Const(e)
    is String -> Var(e)
    is Expression -> e
    else -> error("$e is not an expression")
}

val variables = mutableSetOf<String>()

sealed class Expression {
    infix operator fun plus(other: Any): Add = Add(this, expression(other))
    infix operator fun minus(other: Any): Subtract = Subtract(this, expression(other))
    infix operator fun times(other: Any): Multiply = Multiply(this, expression(other))
    infix operator fun div(other: Any): Divide = Divide(this, expression(other))
    infix fun pow(other: Double): Pow = Pow(this, other)
}

infix operator fun Double.plus(other: Expression): Add = Add(const(this), other)
infix operator fun Double.minus(other: Expression): Subtract = Subtract(const(this), other)
infix operator fun Double.times(other: Expression): Multiply = Multiply(const(this), other)
infix operator fun Double.div(other: Expression): Divide = Divide(const(this), other)

sealed class Binary(
    open val left: Expression,
    open val right: Expression,
    val op: (Double, Double) -> Double,
    private val c: Char
) : Expression() {
    override fun toString(): String = "$left $c $right"
}

class Add(override val left: Expression, override val right: Expression) : Binary(left, right, Double::plus, '+')
class Subtract(override val left: Expression, override val right: Expression) : Binary(left, right, Double::minus, '-')
class Multiply(override val left: Expression, override val right: Expression) : Binary(left, right, Double::times, '*')
class Divide(override val left: Expression, override val right: Expression) : Binary(left, right, Double::div, '/')

data class Pow(val e: Expression, val n: Double) : Expression() {
    override fun toString(): String = "$e^$n"
}

data class Const(val n: Double) : Expression() {
    override fun toString(): String = n.toString()
}

data class Var(val s: String) : Expression() {
    init {
        variables.add(s)
    }

    override fun toString(): String = s
}

fun eval(e: Expression, v: Map<String, Double>): Double = when(e) {
    is Binary -> e.op(eval(e.left, v), eval(e.right, v))
    is Pow -> eval(e.e, v).pow(e.n)
    is Var -> v[e.s] ?: error("no variable in the map")
    is Const -> e.n
}

fun substitute(e: Expression, v: Map<String, Expression>): Expression = when (e) {
    is Add -> substitute(e.left, v) + substitute(e.right, v)
    is Subtract -> substitute(e.left, v) - substitute(e.right, v)
    is Multiply -> substitute(e.left, v) * substitute(e.right, v)
    is Divide -> substitute(e.left, v) / substitute(e.right, v)
    is Pow -> substitute(e.e, v) pow e.n
    is Var -> v[e.s] ?: error("no variable in the map")
    is Const -> e
}

fun contains(e: Expression, v: String): Boolean = when (e) {
    is Binary -> contains(e.left, v) or contains(e.right, v)
    is Pow -> contains(e.e, v)
    is Var -> e.s == v
    is Const -> false
}

fun d(
    e: Binary,
    v: String,
    op: (Expression, Expression) -> Expression
): Expression = if (contains(e, v)) {
    val l = if (contains(e.left, v)) derivative(e.left, v) else e.left
    val r = if (contains(e.right, v)) derivative(e.right, v) else e.right
    op(l, r)
} else const(0.0)

fun derivative(e: Expression, v: String): Expression = when (e) {
    is Add -> derivative(e.left, v) + derivative(e.right, v)
    is Subtract -> derivative(e.left, v) - derivative(e.right, v)
    is Multiply -> d(e, v, Expression::times)
    is Divide -> d(e, v, Expression::div)
    is Var -> if (e.s == v) const(1.0) else const(0.0)
    is Const -> const(0.0)
    is Pow -> if (contains(e.e, v)) {
        val coefficient = if (e.n > 2) e.n - 1 else null
        if (coefficient == null) {
            e.n * e.e
        } else e.n * Pow(e.e, coefficient)
    } else const(0.0)
}
