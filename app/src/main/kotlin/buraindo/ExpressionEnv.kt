package buraindo

import Env
import EvaluationEnv

class ExpressionEnv<TOKEN>(private val evalEnv: EvaluationEnv<TOKEN>) : Env<TOKEN, Expression> {

    override val TOKEN.c: Expression
        get() = const(with(evalEnv) {
            this@c.c
        })

    override val Double.c: Expression
        get() = Const(this)

    override val variableIdentifiers: List<TOKEN>
        get() = evalEnv.variableIdentifiers

    override fun Expression.times(x: Expression): Expression = this.times(x)

    override fun Expression.plus(x: Expression): Expression = this.plus(x)

    override fun Expression.minus(x: Expression): Expression = this.minus(x)

    override fun Expression.div(x: Expression): Expression = this.div(x)

    override fun Expression.pow(x: Double): Expression = this.pow(x)

    override fun Expression.pow(x: Int): Expression = this.pow(x.toDouble())

}