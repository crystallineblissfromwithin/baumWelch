import org.apache.commons.math4.legacy.linear.*;
import org.apache.commons.math4.legacy.optim.MaxEval;
import org.apache.commons.math4.legacy.optim.MaxIter;
import org.apache.commons.math4.legacy.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math4.legacy.optim.univariate.*;
import org.apache.commons.math4.legacy.analysis.*;
import java.util.StringTokenizer;
import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Stream;
/*
Edwin Xie
This program takes in raw data from [ticker].txt, feeds it to Algorithm, then calculates the net loss/gain
 */
public class Trader {
    int numObv;// number of observations the trader makes before it starts trading each day
    int startDate;
    int endDate;
    int N;
    int K;
    int trainT;
    int testT;
    int numoppt = 0;
    int numtrad = 0;
    String ticker;
    ArrayList<ArrayList<Double>> traininglist;
    ArrayList<ArrayList<Double>> testinglist;
    ArrayList<Double> outcomes;
    public Trader()
    {

    }
    public Trader(int traindays,int testdays, int obv, String tick, int n, int k)// dates are YYYYMMDD format
    {
        trainT = traindays;
        testT = testdays;
        ticker = tick;
        numObv = obv;
        N = n;
        K = k;

    }
    public void run()
    {
        try {
            read();
        }
        catch (Exception e)
        {
            System.out.println("reading not done :(");
        }

        Tester test = new Tester(N,K,traininglist,testinglist,outcomes);
        test.test();
        double totalgrowth = 1;
        for(int i = 0;i<test.expGrowth.length;i++)
        {
            totalgrowth*=test.actGrowth[i];
        }
        System.out.println(totalgrowth);
    }
    public void read() throws Exception
    {// makes this class trade the stock 'ticker' by training for trainT, then testing for testT days after that
        BufferedReader br = new BufferedReader(new FileReader("src/^"+ticker+".txt"));
        br.readLine();// throws away first line
        String line = "";
        String[] vals=null;
        int daycount = 0;
        vals = br.readLine().split(",");
        String currentDay=vals[2];
        traininglist = new ArrayList<ArrayList<Double>>();
        while(daycount < trainT)
        {
            traininglist.add(new ArrayList<Double>());
            while (vals[2].equals(currentDay))
            {
                traininglist.get(daycount).add(Double.parseDouble(vals[4]));
                vals = br.readLine().split(",");
            }
            //new day
                daycount++;
                currentDay = vals[2];
        }
        ArrayList<ArrayList<Double>> templist = new ArrayList<ArrayList<Double>>();
        daycount=0;
        while(daycount<testT)
        {
            templist.add(new ArrayList<Double>());
            while (vals[2].equals(currentDay))
            {
                templist.get(daycount).add(Double.parseDouble(vals[4]));
                vals = br.readLine().split(",");
            }
            //new day
            daycount++;
            currentDay = vals[2];
        }
        testinglist = new ArrayList<ArrayList<Double>>();
        outcomes = new ArrayList<Double>();
        for(int i = 0;i<templist.size();i++)
        {
            for(int j = numObv;j<templist.get(i).size()-1;j++)
            {
                ArrayList<Double> temp = new ArrayList<Double>();
                for (int k = 0; k < j; k++)//j represents the length of the list
                {
                    temp.add(templist.get(i).get(k));
                }
                testinglist.add(temp);
                outcomes.add(templist.get(i).get(j));
            }
        }
    }
}
