/*
Edwin Xie
Where the program is run: As of now the growth multipliers for each index are printed out
The data was obtained from stooq.com
Unfortunately this program doesn't consistently generate profits, likely due to the fact that my data is extremely limited
*/
public class Main {
    public static void main(String[] args) {
        String[] tickers = {"aex","aor","ath","bel20","bet","bux","bvp","cac","cdax","dax","djc","dji","djt","dju","fmib","ftm","hex","hsi","ibex","icex","ipc","ipsa","jci","klci","kospi","mdax","moex","mrv","mt30","ndq","ndx","nkx","nomuc","nz50","omxr","omxs","omxt","omxv","oseax","psei","psi20","px","rts","sdxp","set","shbs","shc","smi","sofix","spx","sti","tasi","tdxp","tsx","twse","ukx","xu100"};

        try {
            Trader[] array = new Trader[tickers.length];
            for(int i = 0;i<tickers.length;i++)
            {
                System.out.println("Trading "+tickers[i]+"...");
                array[i] = new Trader(25,7,10,tickers[i],4,5);
                array[i].run();
            }
        }
        catch (Exception e){
            System.out.println("something went wrong!");
        }
    }
}