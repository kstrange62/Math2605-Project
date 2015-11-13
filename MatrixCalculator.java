/**
 * Class to represent Matrix Calculator.
 * @author Kalani Strange
 * @version 1.0
 */
public class MatrixCalculator {
    
    public static float[][] matrixMultiplication(float[][] A, float[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (nA != mB) {
            throw new IllegalArgumentException("Cannot multiply these matrices.");
        }
        float[][] C = new float[mA][nB];
        for (int i = 0; i < mA; i++) {
            for (int j = 0; j < nB; j++) {
                for (int k = 0; k < nA; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }

    public static float[] matrixVectorMultiplication(float[][] A, float[] v) {
        int mA = A.length;
        int nA = A[0].length;
        int nV = v.length;
        if (mA != nV) {
            throw new IllegalArgumentException("Cannot multiply.");
        }
        float[] x = new float[nV];
        for (int i =0; i < mA; i++) {
            for (int j = 0; j < nV; j++) {
                x[i] += A[i][j] * v[j];
            }
        }
        return x;
    }

    public static float[][] scalarMatrixMultiplication(float[][] A, float c) {
        int m = A.length;
        int n = A[0].length;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = A[i][j] * c;
            }
        }
        return A;
    }

    public static float[] scalarVectorMultiplication(float[] v, float c) {
        int n = v.length;
        for (int i = 0; i < n; i++) {
            v[i] = v[i] * c;
        }
        return v;
    }

    public static float[][] addMatrices(float[][] A, float[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (mA != mB || nA != nB) {
            throw new IllegalArgumentException("Cannot add.");
        }
        float[][] C = new float[mA][nA];
        for (int i = 0; i < mA; i++) {
            for (int j = 0; j < nB; j++) {
                C[i][j] = A[i][j] + B[i][j];
            }
        }
        return C;
    }

    public static float[][] subtractMatrices(float[][] A, float[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (mA != mB || nA != nB) {
            throw new IllegalArgumentException("Cannot subtract.");
        }
        float[][] C = new float[mA][nA];
        for (int i = 0; i < mA; i++) {
            for (int j = 0; j < nB; j++) {
                C[i][j] = A[i][j] - B[i][j];
            }
        }
        return C;
    }

    public static float[] subtractVectors(float[] v, float[] x) {
        int nV = v.length;
        int nX = x.length;
        if (nV != nX) {
            throw new IllegalArgumentException("Cannot subtract.");
        }
        float[] u = new float[nV];
        for (int i = 0; i < nV; i++) {
            u[i] = v[i] - x[i];
        }
        return u;
    }

    public static float vectorDotProduct(float[] x, float[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Vectors are not the same length.");
        }
        float sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i] * y[i];
        }
        return sum;
    }

    public static float[][] transpose(float[][] A) {
        int m = A.length;
        int n = A[0].length;
        float[][] C = new float[n][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[j][i] = A[i][j];
            }
        }
        return C;
    }

    public static float[] normalizeVector(float[] v) {
        float[] tmp = new float[v.length];
        float magnitude = 0;
        for (int i = 0; i < v.length; i++) {
            magnitude += v[i] * v[i];
        }
        magnitude = (float) Math.sqrt(magnitude);
        for (int i = 0; i < v.length; i++) {
            tmp[i] = v[i] / magnitude;
        }
        return tmp;
    }

    public static float[][] invert2x2Matrix(float[][] A) {
        float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        float[][] inverse = new float[2][2];
        if (det != 0) {
            float tmp = 1 / det;
            inverse = new float[A.length][A[0].length];
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                    inverse[i][j] = A[i][j] * tmp;
                }
            }
            float tmp0 = inverse[0][0];
            inverse[0][0] = inverse[1][1];
            inverse[1][1] = tmp0;
            inverse[0][1] = inverse[0][1] * -1;
            inverse[1][0] = inverse[1][0] * -1;
        }
        return inverse;
    }

    public static float calculateMatrixTrace(float[][] A) {
        int m = A.length;
        int n = A[0].length;
        if (m != n) {
            throw new IllegalArgumentException("Not a square matrix.");
        }
        float trace = 0;
        int j = 0;
        for (int i = 0; i < A.length; i++) {
            trace += A[i][j];
            j++;
        }
        return trace;
    }
    
    public static int power_method(float[][] matrix, float[] vector, float tolerance, int N) {
        if (matrix.length != matrix[0].length) {
            throw new IllegalArgumentException("The matrix must be nxn.");
        }
        if (vector.length != matrix.length) {
            throw new IllegalArgumentException("The initial guess vector must be of length n.");
        }
        int iterations = 0;
        float[] w = new float[matrix.length];
        w[0] = 1;
        for (int i = 1; i < vector.length; i++) {
            w[i] = 0;
        }
        float tmpTolerance =  Float.MAX_VALUE;
        float[] u = vector;
        float eValue = 0;
        float[] unitEVector = new float[2];
        while (N > 0 && tmpTolerance > tolerance) {
//            System.out.println("Iteration: " + (iterations + 1));
//            System.out.println("u vector: " + Arrays.toString(u));

            //calculates eigenvector
            float[] eVector = matrixVectorMultiplication(matrix, u);

            //normalizes eigenvector
            unitEVector = normalizeVector(eVector);
//            System.out.println("Eigenvector: " + Arrays.toString(eVector));
//            System.out.println("Normalized Eigenvector: " + Arrays.toString(unitEVector));

            //calculates eigenvalue
            float oldEValue = eValue;
            eValue = vectorDotProduct(w, eVector);
            float tmp = 0;
            tmp = vectorDotProduct(w, u);
            if (tmp == 0) {
                tmp = 1;
            }
            eValue = eValue / tmp;
            eValue = Math.abs(eValue);
//            System.out.println("Eigenvalue: " + eValue);

            //compare tolerance
            tmpTolerance = Math.abs((eValue - oldEValue) / eValue);
//            System.out.println("Tolerance: " + tmpTolerance);

            u = eVector;
            iterations++;
            N--;
//            System.out.println();
        }

        //return number of iterations
        if (iterations + 1 >= 99 && tmpTolerance > tolerance) {
            return -1;
        }
        System.out.println("Eigenvalue: " + eValue);
        System.out.println("Eigenvector: " + Arrays.toString(unitEVector));
        return iterations + 1;
    }

    public static ArrayList<float[][]> generate1000_2x2Matrices() {
        ArrayList<float[][]> matrices = new ArrayList<>();
        Random rand = new Random();
        float max = 2;
        for (int  k = 0; k < 1000; k++) {
            float[][] tmp = new float[2][2];
            for (int i = 0; i < tmp.length; i++) {
                for (int j = 0; j < tmp[0].length; j++) {
                    tmp[i][j] = rand.nextFloat() * ((max * -1) - max) + max;
                }
            }
            matrices.add(tmp);
        }
        return matrices;
    }

    public static void main(String[] args) {
        ArrayList<float[][]> matrices = generate1000_2x2Matrices();

        //prints matrices
        System.out.println("Matrices:");
        for (int i = 0; i < 1000; i++) {
            System.out.println("Matrix " + i);
            float[][] matrix = matrices.get(i);
            System.out.println(matrix[0][0] + " " + matrix[0][1]);
            System.out.println(matrix[1][0] + " " + matrix[1][1]);
            System.out.println();
        }

        //check inverse
        ArrayList<float[][]> inverses = new ArrayList<>();
        ArrayList<Float> dets = new ArrayList<>();
        for (float[][] cur: matrices) {
            float det = cur[0][0] * cur[1][1] - cur[0][1] * cur[1][0];
            boolean done = false;
            while (!done) {
                if (det != 0) {
                    float[][] inverse = invert2x2Matrix(cur);
                    inverses.add(inverse);
                    dets.add(det);
                    done = true;
                } else {
                    float max = 2;
                    Random rand = new Random();
                    for (int i = 0; i < cur.length; i++) {
                        for (int j = 0; j < cur[0].length; j++) {
                            cur[i][j] = rand.nextFloat() * ((max * -1) - max) + max;
                        }
                    }
                }
            }
        }

        //prints inverse matrices
        System.out.println("Inverses:");
        for (int i = 0; i < 1000; i++) {
            System.out.println("Inverse of matrix " + i);
            float[][] matrix = inverses.get(i);
            System.out.println(matrix[0][0] + " " + matrix[0][1]);
            System.out.println(matrix[1][0] + " " + matrix[1][1]);
            System.out.println();
        }

        //prints determinants of matrices
        System.out.println("Determinants: ");
        int i = 0;
        for (float cur: dets) {
            System.out.println("Determinant of matrix " + (i++));
            System.out.println(cur);
            System.out.println();
        }

        //power method on A
        float[] u = new float[2];
        u[0] = 1;
        u[1] = 0;
        i = 0;
        for (float[][] cur: matrices) {
            System.out.println("Matrix " + (i++));
            System.out.println("Iterations: " + power_method(cur, u, (float) .00005, 100));
            System.out.println();
        }

        //power method on A inverse
        i = 0;
        for (float[][] cur: inverses) {
            System.out.println("Inverse of matrix " + (i++));
            System.out.println("Iterations: " + power_method(cur, u, (float) .00005, 100));
            System.out.println();
        }

        //traces
        ArrayList<Float> traces = new ArrayList<>();
        float trace = 0;
        for (float[][] cur: matrices) {
            trace = calculateMatrixTrace(cur);
            traces.add(trace);
        }

        //prints traces of matrices
        System.out.println("Traces:");
        for (int k = 0; k < 1000; k++) {
            System.out.println("Trace of matrix " + k);
            System.out.println(traces.get(k));
            System.out.println();
        }
    }
}
