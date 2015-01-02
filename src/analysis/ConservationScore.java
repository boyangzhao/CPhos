package analysis;

import java.util.Iterator;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.awt.Color;

import Jama.Matrix;

/**
 * @author Boyang Zhao
 * @version 1.0
 * @description implementation of algorithms for conservation calculation
 * IC value based entropy calculations
 * results of the analysis are saved in Analysis.cscoresDB
 */

 /**
  * about initialization
  * all the static variables are settings that has to be initialized separately
  * these setting static variables are blank even after instantization of class
  * in this program, the values for these static variables are filled in at start of program with the Utilities.updateSettings() method
  */

public class ConservationScore {

    //internal book keepings
    private char[] aa;
    private class AAGroupings {

        public ArrayList<Character> groupingsLabel;
        public ArrayList<String> groupings;
        public ArrayList<Double> bgFreq;
        public HashMap groupingsAA;

        public AAGroupings(ArrayList<Character> groupingsLabel, ArrayList<String> groupings) {
            this.groupingsLabel = groupingsLabel;
            this.groupings = groupings;
            bgFreq = new ArrayList<Double>();
            groupingsAA = new HashMap<Character, Double>();
            for (int i = 0; i < groupings.size(); i++) {
                Double freqSum = new Double(0);
                for (int j = 0; j < groupings.get(i).length(); j++) {
                    groupingsAA.put(groupings.get(i).charAt(j), i);
                    freqSum += (Double) aaBgFreq.get(groupings.get(i).charAt(j));
                }
                bgFreq.add(freqSum);
            }
        }
    }

    //settings
        //internal
    private HashMap aaBgFreq;
    private boolean bgFreqCalc;
    private double maxEntropy;
    private Analysis analysis;
        //public
    public static double logbase;
    public static ArrayList<Double> bgfreq_a;

    //for sequence logo
    public double[][] icmatrix;
    public double[] icsum;
    public HashMap<Character, Color> groupingsColor;
    public ArrayList<Character> groupingsLabel;
    public ArrayList<String> groupings;
    
    /**
     * construtor for calculating information content/relative entropy only
     * @param bgFreqCalcI
     */
    public ConservationScore(char[][] seqAlign, boolean bgFreqCalcI){
        initialize(bgFreqCalcI);
        //TODO to be finished below is temporary
        /*temporary*/
        //groupingsLabel = new ArrayList(Arrays.asList('K','D','S','N','F','I','G','P','C'));
        //groupings = new ArrayList(Arrays.asList("KRH","DE","ST","NQ","FWY","ILMVA","G","P","C"));
        groupingsLabel = new ArrayList(Arrays.asList('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'));
        groupings = new ArrayList(Arrays.asList("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"));
        groupingsColor = new HashMap<Character,Color>();
        groupingsColor.put('K',Color.BLUE);
        groupingsColor.put('D',Color.CYAN);
        groupingsColor.put('S',Color.GRAY);
        groupingsColor.put('N',Color.MAGENTA);
        groupingsColor.put('F',Color.RED);
        groupingsColor.put('I',Color.YELLOW);
        groupingsColor.put('G',Color.ORANGE);
        groupingsColor.put('P',Color.LIGHT_GRAY);
        groupingsColor.put('C',Color.PINK);
        /*end of temp*/

        AAGroupings AAg = new AAGroupings(groupingsLabel, groupings);
        calculateIC ic = this.new calculateIC(seqAlign, AAg);
        icmatrix = ic.ICMatrix.getArray();
        icsum = ic.ICSum;
    }

    /**
     * constructor for calculating conservation scores from sequence alignment
     * @param analysis
     * @param bgFreqCalcI
     */
    public ConservationScore(Analysis analysisI, boolean bgFreqCalcI){
        analysis = analysisI;
        initialize(bgFreqCalcI);
        calculateCScores();
    }

    private void set_aaBgFreq(){
        aaBgFreq = new HashMap<Character, Integer>();
        for(int i = 0; i < bgfreq_a.size(); i++){
            aaBgFreq.put(aa[i],bgfreq_a.get(i));
        }
    }

    private void initialize(boolean bgFreqCalcI){
        maxEntropy = Utilities.log(20, logbase);
        bgFreqCalc = bgFreqCalcI;
        char[] tmp = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
        aa = tmp;
        set_aaBgFreq();
    }

    /**
     * @modifies Analysis.cscoresDB
     */
    private void calculateCScores(){
        Iterator<char[][]> iter = analysis.phosphositeAlignDB.iterator();
        while(iter.hasNext()){
            int i = analysis.cscoresDB.size();
            ArrayList<String> entry = new ArrayList<String>();
            entry.add(analysis.phosphositeDB.get(i).get(0));
            entry.add(analysis.phosphositeDB.get(i).get(2));
            analysis.cscoresDB.add(entry);
            calculateCScore_single(iter.next());
        }
    }

    private void calculateCScore_single(final char[][] siteAlign){
        ArrayList<Character> groupingsLabel;
        ArrayList<String> groupings;
        AAGroupings AAg1, AAg2, AAg3;
        String ICcomponents;

        //calculate1 {single's}
        groupingsLabel = new ArrayList(Arrays.asList('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'));
        groupings = new ArrayList(Arrays.asList("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"));
        AAg1 = new AAGroupings(groupingsLabel, groupings);
        calculateIC ic1 = this.new calculateIC(siteAlign, AAg1);

        //calculate2 {KRH, DE, ST, NQ, FWY, ILMVA, G, P, C}
        groupingsLabel = new ArrayList(Arrays.asList('K','D','S','N','F','I','G','P','C'));
        groupings = new ArrayList(Arrays.asList("KRH","DE","ST","NQ","FWY","ILMVA","G","P","C"));
        AAg2 = new AAGroupings(groupingsLabel, groupings);
        calculateIC ic2 = this.new calculateIC(siteAlign, AAg2);

        //calculate3 {KRH, STDE, NQ, FWY, ILMVA, G, P, C}
        groupingsLabel = new ArrayList(Arrays.asList('K','S','N','F','I','G','P','C'));
        groupings = new ArrayList(Arrays.asList("KRH","STDE","NQ","FWY","ILMVA","G","P","C"));
        AAg3 = new AAGroupings(groupingsLabel, groupings);
        calculateIC ic3 = this.new calculateIC(siteAlign, AAg3);

        //phosphosite score = sigma(ICs of site only)/(maxEntropy*total num of IC calcs)
        //using ic1, ic2, and ic3
        int sitepos = ic1.ICSum.length/2;
        Double siteScore = (ic1.ICSum[sitepos]+ic2.ICSum[sitepos]+ic3.ICSum[sitepos])/(Utilities.log(AAg1.groupingsLabel.size(), logbase)+Utilities.log(AAg2.groupingsLabel.size(), logbase)+Utilities.log(AAg3.groupingsLabel.size(), logbase));
        analysis.cscoresDB.get(analysis.cscoresDB.size()-1).add(Double.toString(siteScore));
        
        ICcomponents = "[CS:"+Double.toString(Math.round(ic1.ICSum[sitepos]*1e3)/1e3)+"|"+Double.toString(Math.round(ic2.ICSum[sitepos]*1e3)/1e3)+"|"+Double.toString(Math.round(ic3.ICSum[sitepos]*1e3)/1e3)+"]";
        
        //motif score = sigma(IC of all neighboring alignments, excluding site)/(maxEntropy*total num of IC calcs)
        //using ic1 and ic2
        Double motifScore = 0.0;
        int c = 0;
        ICcomponents += "[MS:";
        for(int i = 0; i < ic1.ICSum.length; i++){
            if(i == sitepos || ic1.ICSum[i] == 0.0)
                continue;
            motifScore += ic1.ICSum[i];
            ICcomponents += Double.toString(Math.round(ic1.ICSum[i]*1e3)/1e3)+";";
            c++;
        }
        int c2 = 0;
        ICcomponents = ICcomponents.substring(0, ICcomponents.length()-1);
        ICcomponents += "|";
        for(int i = 0; i < ic2.ICSum.length; i++){
            if(i == sitepos || ic2.ICSum[i] == 0.0)
                continue;
            motifScore += ic2.ICSum[i];
            ICcomponents += Double.toString(Math.round(ic2.ICSum[i]*1e3)/1e3)+";";
            c2++;
        }
        motifScore = motifScore/(Utilities.log(AAg1.groupingsLabel.size(), logbase)*c+Utilities.log(AAg2.groupingsLabel.size(), logbase)*c2);
        ICcomponents = ICcomponents.substring(0, ICcomponents.length()-1);
        ICcomponents += "]";

        analysis.cscoresDB.get(analysis.cscoresDB.size()-1).add(Double.toString(motifScore));
        analysis.cscoresDB.get(analysis.cscoresDB.size()-1).add(ICcomponents);
    }

    /**
     * generates information content score
     */
    private class calculateIC {

        protected double[] ICSum;
        protected Matrix ICMatrix;
        
        public calculateIC(final char[][] siteAlign, final AAGroupings groupings){
            //p(b)
            double[][] freqMatrix_a = new double[siteAlign.length][groupings.groupingsLabel.size()];

            for(int i = 0; i < freqMatrix_a.length; i++){
                for(int j = 0; j < siteAlign[0].length; j++){
                    for(int a = 0; a < aa.length; a++){
                        if(siteAlign[i][j] == aa[a]){
                            freqMatrix_a[i][(Integer) groupings.groupingsAA.get(aa[a])]++;
                            break;
                        }
                    }
                }
            }

            Matrix freqMatrix = new Matrix(freqMatrix_a);
            Matrix probabilityMatrix = freqMatrix.times(1.0/siteAlign[0].length);

            //this is on the side: sum the p(b), if the sum is 0, this is a way to know if everything in the column is blank
            double[] pSum = new double[probabilityMatrix.getRowDimension()];
            for(int i = 0; i < probabilityMatrix.getRowDimension(); i++){
                for(int j = 0; j < probabilityMatrix.getColumnDimension(); j++){
                    pSum[i] += probabilityMatrix.get(i, j);
                }
            }

            //take into account background frequency if applies
            Matrix divided;

            if(bgFreqCalc){
                //p_ref(b)
                double[][] baseFreqMatrix_a2 = new double[siteAlign.length][groupings.groupingsLabel.size()];
                for(int i = 0; i < baseFreqMatrix_a2.length; i++){
                    for(int j = 0; j < baseFreqMatrix_a2[0].length; j++)
                        baseFreqMatrix_a2[i][j] = 1.0/groupings.bgFreq.get(j);
                }
                Matrix baseFreqMatrix = new Matrix(baseFreqMatrix_a2);

                //p(b)/p_ref(b)
                divided = probabilityMatrix.arrayTimes(baseFreqMatrix);
            } else
                divided = probabilityMatrix.copy(); //not really divided, just keep the same variable name for later steps

            //log(p(b)/p_ref(b))
            double[][] lned_a = new double[divided.getRowDimension()][divided.getColumnDimension()];
            for(int i = 0; i < divided.getRowDimension(); i++){
                for(int j = 0; j < divided.getColumnDimension(); j++){
                    double val = divided.get(i, j);
                    lned_a[i][j] = (val == 0) ? 0 : Utilities.log(val, logbase);
                }
            }

            //p(b)*log(p(b)/p_ref(b))
            ICMatrix = probabilityMatrix.arrayTimes(new Matrix(lned_a));

            //sigma(p(b)*log(p(b)/p_ref(b))) -> this gives the relative entropy or information content depending on value of bgFreqCalc
            ICSum = new double[ICMatrix.getRowDimension()];
            for(int i = 0; i < ICMatrix.getRowDimension(); i++){
                for(int j = 0; j < ICMatrix.getColumnDimension(); j++){
                    ICSum[i] += ICMatrix.get(i, j);
                }
            }

            //note: if background ref was not calculated (i.e. bgFreqCalc is false)
            //then the sigma represents the negative of entropy
            //the information content is max entropy (log of 20 base 2) + sigma
            //20 because there are 20 amino acids
            //if pSum is 0.0 for a column, then everything in that column is blank, so IC should be zero, not max entropy
            if(!bgFreqCalc){
                for(int i=0; i < ICSum.length; i++){
                    ICSum[i] = (pSum[i] == 0.0) ? 0.0 : Utilities.log(groupings.groupingsLabel.size(), logbase)+ICSum[i];
                }
            }
        }
    }
}
