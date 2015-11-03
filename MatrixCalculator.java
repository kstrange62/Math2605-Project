/**
 * Class to represent Matrix Calculator.
 * @author Kalani Strange
 * @version 1.0
 */
public class MatrixCalculator {
    public static double[][] matrixMultiplication(double[][] A, double[][] B) {
           int mA = A.length;
           int nA = A[0].length;
           int mB = B.length;
           int nB = B[0].length;
           if (nA != mB) {
               throw new IllegalArgumentException("Cannot multiply these matrices.");
           }
           double[][] C = new double[mA][nB];
           for (int i = 0; i < mA; i++) {
               for (int j = 0; j < nB; j++) {
                   for (int k = 0; k < nA; k++) {
                       C[i][j] += A[i][k] * B[k][j];
                   }
               }
           }
           return C;
        }
        public static double[][] scalarMatrixMultiplication(double[][] A, double c) {
            int m = A.length;
            int n = A[0].length;
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    A[i][j] = A[i][j] * c;
                }
            }
            return A;
        }
        public static double[][] addMatrices(double[][] A, double[][] B) {
            int mA = A.length;
            int nA = A[0].length;
            int mB = B.length;
            int nB = B[0].length;
            if (mA != mB || nA != nB) {
                throw new IllegalArgumentException("Cannot add.");
            }
            double[][] C = new double[mA][nA];
            for (int i = 0; i < mA; i++) {
                for (int j = 0; j < nB; j++) {
                    C[i][j] = A[i][j] + B[i][j];
                }
            }
            return C;
        }
        public static double[][] subtractMatrices(double[][] A, double[][] B) {
            int mA = A.length;
            int nA = A[0].length;
            int mB = B.length;
            int nB = B[0].length;
            if (mA != mB || nA != nB) {
                throw new IllegalArgumentException("Cannot add.");
            }
            double[][] C = new double[mA][nA];
            for (int i = 0; i < mA; i++) {
                for (int j = 0; j < nB; j++) {
                    C[i][j] = A[i][j] - B[i][j];
                }
            }
            return C;
        }
        public static double vectorDotProduct(double[] x, double[] y) {
            if (x.length != y.length) {
                throw new IllegalArgumentException("Vectors are not the same lenght.");
            }
            double sum = 0;
            for (int i = 0; i < x.length; i++) {
                sum += x[i] * y[i];
            }
            return sum;
        }
        public static double[][] transpose(double[][] A) {
            int m = A.length;
            int n = A[0].length;
            double[][] C = new double[n][m];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    C[j][i] = A[i][j];
                }
            }
            return C;
        }
}
