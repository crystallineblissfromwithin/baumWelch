import org.apache.commons.math4.legacy.linear.*;
import java.util.Arrays;
import java.lang.Math;
/*
Edwin Xie
Implementation of Baum-Welch Algorithm on a discrete HMM with multiple sequences of observations
- gamma and xi are not used because of scaling/precision issues
- convergence criterion are as follows:
   - the norm of the change of each of the initial condition matrices/vectors is less than the tolerance
   - the change in the sum of each log of P(Y|theta) is less than the tolerance
*/
public class BaumWelch {
    double tolerance = 0.000001;
    int N;// number of hidden states of X: these states will be labeled 0 to N-1
    int K;// number of possible observations of Y: these states will be labeled 0 to K-1

    int L;// number of sequences
    int[][] Y; // each Y[l] is a sequence of observations that the random variable Y took: Y[l][t] is the observation at time t; t ranges 0 to T[l]-1
    int[] T;// T[l] is length of Y[l]
    double[][] A;// transition matrix: A[i][j] is the transition probability i->j
    double frobA;// frobenius norm (normalized to size of A) of the change of A when updated
    double[][] B;// emissions matrix B[i][k] is the probability, given hidden state i, of observation k
    double frobB;// frobenius norm (normalized to size of B) of the change of B when updated
    double[] P;// 'pi' in literature, the initial probability distribution of X
    double frobP;// frobenius norm (normalized to size of P) of the change of P when updated
    double[][][] Malpha;// stores malpha values
    double[][][] Mbeta;// stores mbeta values
    double[][] C;// scaling factors to avoid precision issues
    double logC;

    public BaumWelch(int n, int k, int[][] y) {
        N = n;
        K = k;
        Y = y;
        L = Y.length;
        T = new int[L];
        for(int l = 0;l<L;l++)
        {
            T[l] = Y[l].length;
        }
        initialize();
        randomize();
    }

    public BaumWelch(int n, int k, int[][] y, double[][] a, double[][] b, double[] p) {
        N = n;
        K = k;
        Y = y;
        L = Y.length;
        T = new int[L];
        for(int l = 0;l<L;l++)
        {
            T[l] = Y[l].length;
        }
        initialize();
        A = a;
        B = b;
        P = p;
    }
    public void train()
    {
        int count = 0;
        double templogC = 0;
        //info();
        count++;
        calculate();
        update();
        // calculate logC
        logC = 0;
        for(int l = 0;l<L;l++)
        {
            for (int t = 0; t < T[l]; t++) {
                logC += Math.log(C[l][t]);
            }
        }
        //System.out.println("Round "+count+" of training successful!");
        //info();
        do
        {
            templogC = logC;

            count++;
            calculate();
            update();
            logC = 0;
            for(int l = 0;l<L;l++)
            {
                for (int t = 0; t < T[l]; t++) {
                    logC += Math.log(C[l][t]);
                }
            }
            if(count>1000000)//for reasonable N and K values this doesn't ever happen - however it takes
            {
                info();
                break;
            }
            //info();
        }while(!(frobA<tolerance && frobB<tolerance && frobP < tolerance && Math.abs(templogC-logC)<tolerance));
        calculate();
        //info();
    }
    public void initialize() {
        A = new double[N][N];
        B = new double[N][K];
        P = new double[N];
        C = new double[L][];
        for(int l = 0;l<L;l++)
        {
            C[l] = new double[T[l]];
        }
        Malpha = new double[L][][];
        Mbeta = new double[L][][];
        for(int l = 0;l<L;l++)
        {
            Malpha[l] = new double[N][T[l]];
            Mbeta[l] = new double[N][T[l]];
        }
        frobA=0;
        frobB=0;
        frobP=0;
    }

    public void randomize() {
        double[] tempX = new double[N + 1];
        for (int i = 0; i < N; i++) {
            tempX[0] = 0;
            for (int j = 1; j < N; j++) {
                tempX[j] = Math.random();
            }
            tempX[N] = 1;
            Arrays.sort(tempX);
            for (int j = 0; j < N; j++) {
                A[i][j] = tempX[j + 1] - tempX[j];
            }
        }
        // randomizes values of B
        double[] tempY = new double[K + 1];
        for (int i = 0; i < N; i++) {
            // doing this to ensure that the randomization is uniform
            tempY[0] = 0;
            for (int j = 1; j < K; j++) {
                tempY[j] = Math.random();
            }
            tempY[K] = 1;
            Arrays.sort(tempY);
            for (int j = 0; j < K; j++) {
                B[i][j] = tempY[j + 1] - tempY[j];
            }
        }
        // randomizes values of P
        tempX[0] = 0;
        for (int j = 1; j < N; j++) {
            tempX[j] = Math.random();
        }
        tempX[N] = 1;
        Arrays.sort(tempX);
        for (int j = 0; j < N; j++) {
            P[j] = tempX[j + 1] - tempX[j];
        }
    }

    public double alpha(int l, int i, int t)//calculates the unscaled alpha value
    {
        double sum = 0;
        for(int j = 0; j<N;j++){
            sum+=Malpha[l][j][t-1]*A[j][i];
        }
        return sum*B[i][Y[l][t]];
    }

    public double beta(int l, int i, int t)//calculates the unscaled beta value
    {
        double sum = 0;
        for(int j = 0;j<N;j++)
        {
            sum+=Mbeta[l][j][t+1]*A[i][j]*B[j][Y[l][t+1]];
        }
        return sum;
    }
    public void calculate()// calculate Malpha, Mbeta, C
    {//scale, memoize and store in Malpha, and do the same for C in the process
        double tempsum;
        for (int l = 0;l<L;l++)
        {
            //calculate first column and C[0]
                tempsum = 0;
            for (int i = 0; i < N; i++) {
                Malpha[l][i][0] = P[i] * B[i][Y[l][0]];
                tempsum += Malpha[l][i][0];
            }
            C[l][0] = 1.0 / tempsum;
            //normalize first column
            for (int i = 0; i < N; i++) {
                Malpha[l][i][0] *= C[l][0];
            }
            //For each column from 1 to T-1:
            for (int t = 1; t < T[l]; t++) {
                //calculate t-th column and C[t]
                tempsum = 0;
                for (int i = 0; i < N; i++) {
                    Malpha[l][i][t] = alpha(l,i, t);
                    tempsum += Malpha[l][i][t];
                }
                C[l][t] = 1.0 / tempsum;
                //normalize first column
                for (int i = 0; i < N; i++) {
                    Malpha[l][i][t] *= C[l][t];
                }
            }
            //scale, memoize and store Mbeta
            //calculate last column
            for (int i = 0; i < N; i++) {
                Mbeta[l][i][T[l] - 1] = C[l][T[l] - 1];
            }
            //use column t to calculate column t-1
            for (int t = T[l] - 2; t >= 0; t--) {
                for (int i = 0; i < N; i++) {
                    Mbeta[l][i][t] = beta(l,i, t) * C[l][t];
                }
            }
        }
    }
    public void update()// calculate and update A,B,P
    {
        // calculate denominators of A and B together because they are similar
        double[] denomA = new double[N];
        double[] denomB = new double[N];
        double tempsum;
        for(int i = 0;i<N;i++)
        {
            for (int l=0;l<L;l++)
            {
                tempsum = 0;
                for (int t = 0; t < T[l] - 1; t++) {
                    tempsum += Malpha[l][i][t] * Mbeta[l][i][t] / C[l][t];
                }
                denomA[i]+=tempsum;
                denomB[i] += (tempsum + (Malpha[l][i][T[l] - 1] * Mbeta[l][i][T[l] - 1] / C[l][T[l] - 1]));
            }
        }
        // calculate and update A and frobA
        frobA = 0;
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                tempsum = 0;
                for(int l=0;l<L;l++)
                {
                    for (int t = 0; t < T[l] - 1; t++)
                    {
                        tempsum += (Malpha[l][i][t] * B[j][Y[l][t + 1]] * Mbeta[l][j][t + 1]);
                    }
                }
                frobA += Math.pow(A[i][j]*(1-(tempsum/denomA[i])),2);
                A[i][j]*=(tempsum/denomA[i]);
            }
        }
        frobA = Math.pow(frobA/(Math.pow(N,2)),0.5);
        // calculate and update B
        frobB = 0;
        for(int i = 0;i<N;i++)
        {
            for(int k = 0;k<K;k++)
            {
                tempsum = 0;
                for(int l=0;l<L;l++)
                {
                    for (int t = 0; t < T[l]; t++) {
                        tempsum += kronk(Y[l][t], k) * Malpha[l][i][t] * Mbeta[l][i][t] / C[l][t];
                    }
                }
                frobB+=Math.pow(B[i][k]-(tempsum/denomB[i]),2);
                B[i][k]=tempsum/denomB[i];
            }
        }
        frobB = Math.pow(frobB/(N*K),0.5);
        // calculate and update P
        frobP = 0;

        double temparr[] = new double[L];
        for(int l =0;l<L;l++)
        {
            for (int j = 0; j < N; j++) {
                temparr[l] += Malpha[l][j][0] * Mbeta[l][j][0];
            }
        }
        for(int i =0;i<N;i++)
        {
            tempsum=0;
            for(int l=0;l<L;l++)
            {
                tempsum+= Malpha[l][i][0] * Mbeta[l][i][0]/temparr[l];
            }
            tempsum/=L;
            frobP+=Math.pow(tempsum-P[i],2);
            P[i] = tempsum;
        }
        frobP = Math.pow(frobP/N,0.5);
    }
    public int kronk(int m, int n)// kronecker delta
    {
        return (m == n) ? 1 : 0;
    }
    public void info()//debug method
    {
        for(int i = 0;i<N;i++)
        {
            for(int j = 0;j<N;j++)
            {
                System.out.print(A[i][j]+" ");
            }
            System.out.println();
        }
        System.out.println();
        for(int i = 0;i<N;i++)
        {
            for(int k = 0;k<K;k++)
            {
                System.out.print(B[i][k]+" ");
            }
            System.out.println();
        }
        System.out.println();
        for(int i=0;i<N;i++)
        {
            System.out.print(P[i]+" ");
        }
        System.out.println();
        System.out.println();
        System.out.println(frobA+","+frobB+","+frobP);
    }
}

