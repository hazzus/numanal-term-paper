package buraindo

import Matrix
import System
import Vector
import dot
import inverted
import times
import transposed

fun f(system: System, v: Map<String, Double>): Vector = system.map { eval(it, v) }

// TODO: метод Ньютона
/*fun solve(
    system: System,
    jacobi: Matrix<Expression>,
    initialApproach: Vector
): Sequence<Vector> =
    generateSequence(initialApproach) { x_k ->
        val v = variables.zip(x_k).toMap()
        val fx = f(system, v)
        val matrix = jacobi.map {
            f(it, v)
        }.t().inv()
        x_k.zip(matrix * fx, Double::minus)
    }*/

//
fun solveUniversal(
    F: Expression,
    system: System,
    jacobi: Matrix<Expression>,
    initialApproach: Vector
): Sequence<Vector> =
    generateSequence(initialApproach) { x_k ->
        val v = variables.zip(x_k).toMap()
        val fx = f(system, v)
        val matrix = jacobi.map {
            f(it, v)
        }.transposed().inverted()
        val delta = (matrix * fx).map { varOf("alpha") * it }
        val arg = x_k.zip(delta, Double::minus)
        val q = substitute(F, variables.zip(arg).toMap()) // F(x_k + alpha * delta x_k)
        x_k // TODO: здесь найти минимум q
    }

// TODO: решает систему универсальным методом, раскомментировать, когда метод будет готов
/*fun solveChemicalsUniversal() {
    val t = 623.15
    val res = computeCoefficients(t)
    val k1 = res.K1
    val k2 = res.K2
    val k3 = res.K3
    val dHCl = res.dHCl
    val dH2 = res.dH2
    val dAlCl = res.dAlCl
    val dAlCl2 = res.dAlCl2
    val dAlCl3 = res.dAlCl3
    val pAlCl = varOf("p_alcl") //
    val pH2 = varOf("p_h2") //
    val pHCl = varOf("p_hcl") //
    val pAlCl2 = varOf("p_alcl2") //
    val pAlCl3 = varOf("p_alcl3") //
    val pGAlCl = res.pGAlCl
    val pGHCl = res.pGHCl
    val pGH2 = res.pGH2
    val pGAlCl2 = res.pGAlCl2
    val pgAlCl3 = res.pGAlCl3

    val system = listOf(
        k1 * pAlCl.pow(2.0) * pH2 - pHCl.pow(2.0),
        k2 * pAlCl2 * pH2 - pHCl.pow(2.0),
        k3 * pAlCl3.pow(2.0) * pH2.pow(3.0) - pHCl.pow(6.0),
        dHCl * (pGHCl - pHCl) + 2.0 * dH2 * (pGH2 - pH2),
        dAlCl * (pGAlCl - pAlCl) + 2.0 * dAlCl2 * (pGAlCl2 - pAlCl2) + 3.0 * dAlCl3 * (pgAlCl3 - pAlCl3) + dHCl * (pGHCl - pHCl)
    )

    val jacobi = variables.map { v -> system.map { e -> derivative(e, v) } }
    val quadraticForm = system.foldRight(const(0.0)) { e: Expression, acc: Expression -> (e pow 2.0) + acc }
    val initialApproach = listOf(-0.1, -0.1, -0.1, -0.1, -0.1)
    val eps = 1e-6
    val solution = solveUniversal(quadraticForm, system, jacobi, initialApproach).takeWhile { v ->
        f(system, variables.zip(v).toMap()).run { dot(this, this) } > eps
    }
    solution.last().zip(variables).forEach { (value, name) ->
        println("$name: $value")
    }
    print(res.print())
}*/

// TODO: решает систему методов Ньютона
/*fun solveChemicals() {
    val t = 623.15
    val res = computeCoefficients(t)
    val k1 = res.K1
    val k2 = res.K2
    val k3 = res.K3
    val dHCl = res.dHCl
    val dH2 = res.dH2
    val dAlCl = res.dAlCl
    val dAlCl2 = res.dAlCl2
    val dAlCl3 = res.dAlCl3
    val pAlCl = varOf("p_alcl") //
    val pH2 = varOf("p_h2") //
    val pHCl = varOf("p_hcl") //
    val pAlCl2 = varOf("p_alcl2") //
    val pAlCl3 = varOf("p_alcl3") //
    val pGAlCl = res.pGAlCl
    val pGHCl = res.pGHCl
    val pGH2 = res.pGH2
    val pGAlCl2 = res.pGAlCl2
    val pgAlCl3 = res.pGAlCl3

    val system = listOf(
        k1 * pAlCl.pow(2.0) * pH2 - pHCl.pow(2.0),
        k2 * pAlCl2 * pH2 - pHCl.pow(2.0),
        k3 * pAlCl3.pow(2.0) * pH2.pow(3.0) - pHCl.pow(6.0),
        dHCl * (pGHCl - pHCl) + 2.0 * dH2 * (pGH2 - pH2),
        dAlCl * (pGAlCl - pAlCl) + 2.0 * dAlCl2 * (pGAlCl2 - pAlCl2) + 3.0 * dAlCl3 * (pgAlCl3 - pAlCl3) + dHCl * (pGHCl - pHCl)
    )

    val jacobi = variables.map { v -> system.map { e -> derivative(e, v) } }
    val initialApproach = listOf(-0.1, -0.1, -0.1, -0.1, -0.1)
    val eps = 1e-6
    val solution = solve(system, jacobi, initialApproach).takeWhile { v ->
        f(system, variables.zip(v).toMap()).run { dot(this, this) } > eps
    }
    solution.last().zip(variables).forEach { (value, name) ->
        println("$name: $value")
    }
    print(res.print())
}*/

// TODO: решает простейшую систему универсальным методом (еще не готов)
fun solveTrivial() {
    val x = varOf("x")
    val y = varOf("y")
    val z = varOf("z")
    val system = listOf(
        (x pow 2.0) - 2.0 * x + (y pow 2.0) - z + 1.0,
        x * (y pow 2.0) - x - 3.0 * y + y * z + 2.0,
        x * (z pow 2.0) - 3.0 * z + y * (z pow 2.0) + x * y
    )
    val jacobi = variables.map { v -> system.map { e -> derivative(e, v) } }
    val quadraticForm = system.foldRight(const(0.0)) { e: Expression, acc: Expression -> (e pow 2.0) + acc }
    val initialApproach = listOf(1.0, 1.1, 123.99)
    val eps = 1e-25
    val solution = solveUniversal(quadraticForm, system, jacobi, initialApproach).takeWhile { v ->
        f(system, variables.zip(v).toMap()).run { dot(this, this) } > eps
    }
    solution.last().zip(variables).forEach { (value, name) ->
        println("$name: $value")
    }
}

fun main() {
    solveTrivial()
}