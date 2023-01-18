import org.apache.commons.math4.legacy.analysis.function.Abs;
import org.apache.commons.math4.legacy.linear.*;
import org.apache.commons.math4.legacy.optim.MaxEval;
import org.apache.commons.math4.legacy.optim.MaxIter;
import org.apache.commons.math4.legacy.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math4.legacy.optim.univariate.*;
import org.apache.commons.math4.legacy.analysis.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.function.DoubleUnaryOperator;
/*
Edwin Xie
Given raw data, this class parses it into percentage changes, uses (a generalized version of) the Kelly fraction to find the percentage of the portfolio to "bet" on the stock, and then bets it.
 */
public class Algorithm {
    final double epsilon =1e-14;
    int N;
    int K;
    int[][] Y;
    double m,M;
    boolean manual;
    ArrayList<ArrayList<Double>> raw;
    double[] bins;
    double[][] Delta;
    BaumWelch model;
    public Algorithm(int n, int k, ArrayList<ArrayList<Double>> rawdata, double min, double max)// N and K are the same parameters as in BaumWelch, rawdata is the raw stock prices
    {
        N = n;
        K = k;
        raw = rawdata;
        m = min;
        M = max;
        bins = new double[k+1];
        manual = false;
    }

    public Algorithm(int n, int k, ArrayList<ArrayList<Double>> rawdata)// N and K are the same parameters as in BaumWelch, rawdata is the raw stock prices
    {
        N = n;
        K = k;
        raw = rawdata;
        bins = new double[k+1];
        manual = true;
    }

    public void parse()// parses raw data, creates equally spaced bins, creates Y. Then trains an HMM model on that data.
    {
        Delta = new double[raw.size()][];
        ArrayList<Double> temp = new ArrayList<Double>();
        for(int i =0;i<raw.size();i++)
        {
            Delta[i] = new double[raw.get(i).size()-1];
            for(int j = 0;j<raw.get(i).size()-1;j++)
            {
                Delta[i][j] = raw.get(i).get(j+1)/raw.get(i).get(j);
                //System.out.print(Delta[i][j]+":");
                temp.add(Delta[i][j]);
            }
            //System.out.println();
        }
        Collections.sort(temp);

        if(manual)
        {
            m = Double.POSITIVE_INFINITY;
            M = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < Delta.length; i++) {
                for (int j = 0; j < Delta[i].length; j++) {
                    m = Math.min(m, Delta[i][j]);
                    M = Math.max(M, Delta[i][j]);
                }
            }
        }

        for(int i = 0;i<=K;i++)
        {
            bins[i] = (m*(K-i)+M*(i))/K;
            //System.out.print(bins[i]+":");
        }

        Y = new int[Delta.length][];
        for(int i = 0;i<Delta.length;i++)
        {
            Y[i] = new int[Delta[i].length];
            for(int j = 0;j<Y[i].length;j++)
            {
                for (int k = K; k >= 1; k--) {
                    if (Delta[i][j] < bins[k] + epsilon) {
                        Y[i][j] = k - 1;
                    }
                }
            }
        }
    }
    public void train()
    {
        BaumWelch[] candidates = new BaumWelch[20];
        double[] fitness = new double[candidates.length];// -log(P(Y[0],Y[1]...|theta)), we actually want to minimize this
        double maxfitness = Double.POSITIVE_INFINITY;
        int mostfit = -1;
        for(int i = 0;i<candidates.length;i++) {
            candidates[i] = new BaumWelch(N, K, Y);
            candidates[i].train();
            fitness[i] = candidates[i].logC;
            //System.out.println(i+"fitness:"+fitness[i]);
            if(fitness[i]<maxfitness)
            {
                maxfitness = fitness[i];
                mostfit = i;
            }
            //System.out.println(i+"done!");
        }
        model = new BaumWelch(N,K,Y,candidates[mostfit].A, candidates[mostfit].B,candidates[mostfit].P);
    }
    public UnivariatePointValuePair kelly(int[] O)// given an observation sequence O, decides whether to trade something, returns the fraction of the portfolio to trade on it
    {
        double[] distr = predict(model, O);
        //calculates conservative estimate for EV, returns -1 immediately if EV is negative
        double ev = 0;
        for(int i =0;i<K;i++) {
            ev += ((bins[i]+bins[i+1])/2-1) * distr[i];
        }
        if(ev<0)
        {
            //System.out.println("EV: "+ev);
            return new UnivariatePointValuePair(-1,0);
        }
        EVFunction temp = new EVFunction(distr,bins);
        UnivariateObjectiveFunction func= new UnivariateObjectiveFunction(temp);
        BrentOptimizer bo = new BrentOptimizer(1e-10,1e-14);
        UnivariatePointValuePair max = bo.optimize(GoalType.MAXIMIZE,new SearchInterval(0,1,0.5), func, MaxEval.unlimited(), MaxIter.unlimited());
        return max;
    }
    public static double[] predict(BaumWelch theta, int[] O)
    {
        for(int i = 0;i<O.length;i++)
        {
            //System.out.print(O[i]+";");
        }
        int[][] Y = new int[1][];
        Y[0] = O;
        BaumWelch pred = new BaumWelch(theta.N,theta.K,Y,theta.A,theta.B,theta.P);

        pred.calculate();// now we have the Malpha and Mbeta values we need

        // calculate Gamma_i(T)
        double tempsum = 0;
        for(int i = 0;i<pred.N;i++)
        {
            tempsum+=pred.Malpha[0][i][pred.T[0]-1]*pred.Mbeta[0][i][pred.T[0]-1];
        }
        double[] gammavector = new double[pred.N];
        for(int i = 0;i<pred.N;i++)
        {
            gammavector[i] = pred.Malpha[0][i][pred.T[0]-1]*pred.Mbeta[0][i][pred.T[0]-1]/tempsum;
        }
        RealMatrix transform = MatrixUtils.createRealMatrix(theta.B).transpose().multiply(MatrixUtils.createRealMatrix(theta.A).transpose());
        double[] val = transform.operate(gammavector);
        /*
        for(int i = 0;i<val.length;i++)
        {
           System.out.print(val[i]+":");
        }
        System.out.println();
        */
        return val;
    }
}
class EVFunction implements UnivariateFunction
{
    public
    double[] pDistr;
    double[] expVal;
    public EVFunction(double[] pdistr, double[] eval)
    {
     pDistr = pdistr;
     expVal = eval;
    }

    public double value(double v) {
        double val = 0;
        for(int i = 0;i<pDistr.length;i++)
        {
            val+=pDistr[i]*Math.log(1+v*((expVal[i]+expVal[i+1])/2-1));
        }
        return val;
    }

    @Override
    public double applyAsDouble(double x) {
        return UnivariateFunction.super.applyAsDouble(x);
    }

    @Override
    public DoubleUnaryOperator compose(DoubleUnaryOperator before) {
        return UnivariateFunction.super.compose(before);
    }

    @Override
    public DoubleUnaryOperator andThen(DoubleUnaryOperator after) {
        return UnivariateFunction.super.andThen(after);
    }
}
