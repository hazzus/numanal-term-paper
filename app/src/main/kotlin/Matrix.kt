import buraindo.Expression

typealias Vector = List<Double>
typealias Matrix<T> = List<List<T>>
typealias System = List<Expression>

fun dot(n: Vector, m: Vector): Double = n.zip(m, Double::times).sum()

operator fun Matrix<Double>.times(v: Vector): Vector = map { x -> dot(x, v) }

fun mult(a: Matrix<Double>, b: Matrix<Double>): Matrix<Double> =
    b.transposed().map { a * it }.transposed()

fun Matrix<Double>.inverted(): Matrix<Double> {
    val matrix = map(Vector::toMutableList).toMutableList()
    matrix.forEachIndexed { index, mutableList ->
        mutableList.addAll(
            generateSequence(0.0, { 0.0 }).take(size).toMutableList().also { lst -> lst.set(index, 1.0) })
    }
    val n = size
    val processed = mutableListOf<Int>()
    repeat(size) { _ ->
        val i = (0 until n).find { matrix[it][it] != 0.0 && it !in processed }!!
        processed.add(i)
        val divisionCoefficient = matrix[i][i]
        (0 until (n * 2)).forEach { k ->
            matrix[i][k] /= divisionCoefficient
        }
        (0 until n).forEach { j ->
            if (j != i) {
                val subtractionCoeff = matrix[j][i] / matrix[i][i]
                (0 until n * 2).forEach { k ->
                    matrix[j][k] -= matrix[i][k] * subtractionCoeff
                }
            }
        }
    }
    return matrix.map { list -> list.drop(size) }
}

fun <T> Matrix<T>.transposed(): Matrix<T> {
    val result = mutableListOf<MutableList<T>>()
    repeat(first().size) { result.add(mutableListOf()) }
    forEach {
        it.forEachIndexed { index, elem -> result[index].add(elem) }
    }
    return result
}
