package analysis;

import java.io.*;
import java.util.*;

/**
 * @author Boyang Zhao
 * @version 1.0
 * @description implementation of algorithms for analyzing phosphopeptides
 *
 */

 /**
  * about initialization
  * all the static variables are settings that has to be initialized separately
  * these setting static variables are blank even after instantization of class
  * protFileName and orthoFileName are the only two static variables that are initialized during instantization
  * in this program, the values for these static variables are filled in at start of program with the Utilities.updateSettings() method
  */

 /**
  * locus is a deprecated term. the program was originally designed with swiss-prot database in mind but has since switched to HomoloGene
  * locus below represents the unique identifier for the orthologous group
  * definitions:
  * locus = HID|gene symbols
  * locus_species = HID|gene symbols_[species with spaces]
  */

public class Analysis {
    /* data structure declaration */

    //data structures have same total row # - entries correspond
    //{phosphositeSList, phosphositeDB, phosphositeAlignDB, cscoresDB}
    //{locusSList, locusList}

    //databases
    public ArrayList<String> peptideInputs; //[input seq (w/ symbols)]
    public ArrayList<ArrayList<String>> peptideDB; //[[input seq (w/ symbols), phosphosites]]
    public ArrayList<String> peptideDBambi; //[[input seq (w/ symbols)]]
    public ArrayList<ArrayList<String>> protDB; //[[locus, locus_species, accession number, gene symbol, protein seq]]
    public ArrayList<ArrayList<String>> phosphositeDB; //[[locus, locus_species, phosphosite, gene symbol, corresponding phosphosites of other species]]
    public ArrayList<ArrayList<String>> locusList; //[[locus, # of orthologs, start index in orthologDB, end index in orthologDB]]
    public ArrayList<ArrayList<String>> orthologDB; //[[locus, locus_species, accession number, gene symbol, protein seq]]
    public ArrayList<ArrayList<String>> cscoresDB; //[[locus, phosphosite, site conservation score, motif conservation score]]

    //book keepings
        //public
    public ArrayList<char[][]> phosphositeAlignDB; //[[2D residues alignment]]; index of ArrayList corresponds to phophosite, same ordering as phosphositeDB
                                                   //for char 2D array: each row represents each alignment pos (-1, 0, +1) and each col represents different orthologs
    public long runningTime;
    
        //internal
    private TreeSet<String> locusSList; //[locus]
    public HashMap<String, String> phosphositeSList; //key: [locus \t locus_species \t phosphosite \t gene symbol]; value: input peptide indices (index starts at 1)
    private long startTime;
    private long endTime;

    //settings
        //public
    public String inputSpecies;
    public static boolean display;
        //internal
    private final int siteAlignTotalLen;
    private javax.swing.JFrame parentFrame;
    private boolean db_local;
        //paths
    public static String dbFolderPath;
    public static String protFileName;
    public static String orthoFileName;
    public static String exeFolderPath;
    public static String alnOutFolderPath;
    public static ArrayList<String> db_species;
    public static ArrayList<String> db_target;
        //analysis parameters
    public static int siteAlignNgbrLen;
    public static char cPlaceholder; //TODO currently not used
    public static char nPlaceholder; //TODO currently not used

    //error handling
    public boolean analyzeStatus;

    /**
     * Analysis constructor, initialize all variables
     * @param speciespath
     * @param dbpath
     */
    public Analysis(String speciespath, String dbpath, javax.swing.JFrame parent, boolean displayinput) {
        //inializations
        parentFrame = parent;
        display = displayinput;
        peptideInputs = new ArrayList<String>();
        peptideDB = new ArrayList<ArrayList<String>>();
        peptideDBambi = new ArrayList<String>();
        protDB = new ArrayList<ArrayList<String>>();
        phosphositeDB = new ArrayList<ArrayList<String>>();
        locusList = new ArrayList<ArrayList<String>>();
        orthologDB = new ArrayList<ArrayList<String>>();
        phosphositeAlignDB = new ArrayList<char[][]>();
        cscoresDB = new ArrayList<ArrayList<String>>();
        locusSList = new TreeSet<String>();
        phosphositeSList = new HashMap<String, String>();
        protFileName = speciespath;
        orthoFileName = dbpath;
        db_local = orthoFileName.endsWith(".ws") ? false : true;
        analyzeStatus = true;
        siteAlignTotalLen = siteAlignNgbrLen*2+1;
    }
    
    /**
     * dispatch input into analysis methods for computation
     */
    public void analyzePeptides() {
        
        startTime = System.currentTimeMillis();
        
        if (peptideInputs.isEmpty()) {
            Utilities.displayMessage("Please enter in at least one sequence before continue.", parentFrame);
        } else {
            try {
                Utilities.setOverallMax(100);
                createProtDB();
                Utilities.setOverallProgress(25);
                createOrthologDB();
                Utilities.setOverallProgress(50);
                createSeqAlignment();
                Utilities.setOverallProgress(75);
                createCScoresDB();
                Utilities.setOverallProgress(100);
            } catch (Exception e){
                Utilities.displayMessage(e.getMessage(), parentFrame);
                analyzeStatus = false;
            }
        }
        
        endTime = System.currentTimeMillis();
        runningTime = endTime - startTime;
    }

    private void createPhosphositeDB(){
        for(String psite : phosphositeSList.keySet()){
            System.out.print(psite);
            phosphositeDB.add(new ArrayList(Arrays.asList(psite.split("\\t"))));
        }
    }

    /**
     * create an 2-D array of homologous proteins in same-species
     * @modifies peptideDB, locusSList, phosphositeSList, protDB, peptideDBambi
     */
    private void createProtDB() throws Exception {

        Utilities.startTask("Task 1 of 4: Generating Protein Database - parsing local binary database");

        //create an input stream from local binary database file (same species)
        //then save the database in the form of a 2D array
        //structure: [[locus_species],[accession number],[gene symbol],[sequence]]
        ArrayList<ArrayList<String>> db2DArrayList = Utilities.deserializeDB(dbFolderPath+protFileName, 200000);

        if (analyzeStatus) {
            Utilities.startTask("Task 1 of 4: Generating Protein Database", peptideInputs.size());

            inputSpecies = Utilities.getSpecies(db2DArrayList.get(0).get(0), false);

            //start searching for phosphoprotein (same species)
            for (int i = 0; i < peptideInputs.size(); i++) {
                //note input already has all char converted to upper case
                String rawSeq = peptideInputs.get(i);
                String pOnlySeq = rawSeq.replaceAll("[^A-Z\\*]*", ""); //keep the * for phosphorylation, remove all other symbols
                String cleanedSeq = rawSeq.replaceAll("[^A-Z]*", "");

                //find matches of peptide sequence in protein database of same species
                //protein db is parsed as 2D array and saved in db2DArrayList
                ArrayList<String> firstMatch = new ArrayList<String>();
                int matchN = 0;

                int z=0;

                Iterator<ArrayList<String> > iterL = db2DArrayList.iterator();
                while(iterL.hasNext()){z++;
                    boolean match = false;
                    ArrayList<String> protein = iterL.next();
                    match = protein.get(3).toUpperCase().contains(cleanedSeq);

                    if (match) {
                        if (matchN == 0) {
                            firstMatch.add(rawSeq);
                            firstMatch.add(protein.get(0)); //locus_species
                            firstMatch.add(protein.get(1)); //accession number
                            firstMatch.add(protein.get(2)); //gene symbol
                            firstMatch.add(protein.get(3)); //protein sequence
                        } else {
                            //if multiple matches found, save into separate database
                            //can also save this to database; currently there's no use for saving this.
                            //[a 2D array].add(new ArrayList(Arrays.asList(rawSeq,protein.get(0),protein.get(1),protein.get(2),protein.get(3))));
                        }

                        matchN++;
                    }
                }
                
                //if only one corresponding protein is found, then this will be further analyzed
                if (matchN == 1) {
                    //determine phosphosite(s) from this phosphopeptide
                    //first pos number of input seq
                    int indexPeptSeq = firstMatch.get(4).indexOf(cleanedSeq);
                    int checkasterisk = pOnlySeq.indexOf("*", 0);
                    int findasterisk = checkasterisk;
                    int pN = 0; //keeps track of total number of *, need this to calculate the correct site position
                    String phosphosites = "";
                    if (checkasterisk != -1) {
                        //phosphorylation exists, get locus first
                        //get locus and prot to protDB is haven't done so already
                        String locussp = firstMatch.get(1);
                        String locus = Utilities.getLocus(locussp);

                        if (locusSList.add(locus)) {
                            //[locus, locus_species, accession number, gene symbol, protein seq]
                            protDB.add(new ArrayList(Arrays.asList(locus,firstMatch.get(1),firstMatch.get(2),firstMatch.get(3),firstMatch.get(4))));
                        }

                        //now determine the phosphosite positions
                        do {
                            pN++;
                            int phosphosite = indexPeptSeq + findasterisk - pN;
                            char aaphosphosite = firstMatch.get(4).charAt(phosphosite);
                            int realphosphosite = phosphosite + 1; //index starts at 0, but real count starts at 1
                            String hashkey = locus+"\t"+firstMatch.get(1)+"\t"+String.valueOf(aaphosphosite)+realphosphosite+"\t"+firstMatch.get(3);
                            String hashvalue = phosphositeSList.get(hashkey); //hash value keeps track of input indices (index starts at 1) for given phosphosite
                            if(hashvalue != null){
                                hashvalue = hashvalue+";"+String.valueOf(i+1);
                            } else {
                                hashvalue = String.valueOf(i+1);
                            }
                            phosphositeSList.put(hashkey, hashvalue);
                            phosphosites += String.valueOf(aaphosphosite)+realphosphosite+", ";
                            findasterisk++;

                            findasterisk = pOnlySeq.indexOf("*", findasterisk);
                        } while (findasterisk != -1);

                        phosphosites = phosphosites.substring(0, phosphosites.length()-2);

                        peptideDB.add(new ArrayList(Arrays.asList(rawSeq,phosphosites)));
                    } else {
                        peptideDBambi.add(rawSeq);
                        peptideDB.add(new ArrayList(Arrays.asList(rawSeq,"no phosphosite found - see help for more info")));
                    }

                } else if (matchN == 0) {
                    peptideDBambi.add(rawSeq);
                    peptideDB.add(new ArrayList(Arrays.asList(rawSeq,"no protein match found - see help for more info")));
                } else {
                    peptideDBambi.add(rawSeq);
                    peptideDB.add(new ArrayList(Arrays.asList(rawSeq,"multiple protein matches found - see help for more info")));
                }

                Utilities.setTaskProgress(i+1);
            }
        }
    }

    /**
     * create an 2-D array of orthologs across multiple species
     * @modifies orthologDB
     */
    private void createOrthologDB() throws Exception {
        if(db_local)
            Utilities.startTask("Task 2 of 4: Generating Ortholog Database - parsing binary local database");
        else
            Utilities.startTask("Task 2 of 4: Generating Ortholog Database - parsing homologene database");

        if(locusSList.isEmpty())
            return;

        ArrayList<ArrayList<String>> db2DArrayList = new ArrayList<ArrayList<String>>();
        HashMap<String, HashSet> hDB_hash = new HashMap<String, HashSet>();
        if(db_local){
            //create an input stream from local binary database file (multiple species)
            //then save the database in the form of a 2D array
            //structure: [[locus],[accession number],[gene name],[sequence]]
            db2DArrayList = Utilities.deserializeDB(dbFolderPath+orthoFileName, 500000);
        } else {
            //get the hash table (key=HID; value=list of proteins delimited by "|")
            ObjectInputStream inputStream = null;
            inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(Analysis.dbFolderPath + Utilities.gethDBfilename()), 500000));
            inputStream.readObject();
            hDB_hash = (HashMap<String, HashSet>) inputStream.readObject();
            if (inputStream != null) {
                inputStream.close();
            }
        }
        
        //search for proteins from other species that have the same locus
        if (analyzeStatus) {
            Utilities.startTask("Task 2 of 4: Generating Ortholog Database", locusSList.size());

            Iterator<String> iter = locusSList.iterator();
            int i = 0; //track locus
            int start = 0, end = 0; //track index in orthologDB
            while(iter.hasNext()){
                String locus = iter.next();
                int matchesNum = 0;

                if(db_local){
                    //loop through all proteins in db of other species (db2DArrayList)
                    //and pull out all protein sequences that are of the same locus
                    for (ArrayList<String> protein : db2DArrayList) {
                        boolean match = false;
                        match = protein.get(0).startsWith(locus+"_");

                        if(match){
                            //add to orthologDB
                            ArrayList<String> orthMatch = new ArrayList<String>();
                            orthMatch.add(locus);
                            orthMatch.add(protein.get(0)); //locus_species
                            orthMatch.add(protein.get(1)); //accession number
                            orthMatch.add(protein.get(2)); //gene symbol
                            orthMatch.add(protein.get(3)); //protein sequence

                            if(!protein.get(0).contains(inputSpecies)){
                                orthologDB.add(orthMatch);
                                matchesNum++;
                            }
                        }
                    }
                } else {
                    //using NCBI web service
                    String protlist = Utilities.getFromSet(hDB_hash.get(Utilities.getHID(locus)),false);
                    String[] list = protlist.split("\\|\\|");
                    String[] match = new String[5];

                    String taxonomyclass = Utilities.getDBName(orthoFileName, false);
                    taxonomyclass = taxonomyclass.toLowerCase();
                    for(int m = 0; m < list.length; m++){
                        if(!list[m].isEmpty()){
                            match = Utilities.callws_getgp(list[m], taxonomyclass);
                            if(match[0].equals("found")){
                                //add to orthologDB
                                ArrayList<String> orthMatch = new ArrayList<String>();
                                orthMatch.add(locus);
                                orthMatch.add(locus+"_["+match[1]+"]"); //locus_species
                                orthMatch.add(match[2]); //accession number
                                orthMatch.add(match[3]); //gene symbol
                                orthMatch.add(match[4]); //protein sequence

                                if(!match[1].contains(inputSpecies)){
                                    orthologDB.add(orthMatch);
                                    matchesNum++;
                                }
                            }
                        }
                    }
                }

                /**
                 * add input species using the existing protDB
                 * instead of using the target database (this avoids missing information/version difference
                 * of input between between protDB and target database)
                 */
                for(int p = 0; p < protDB.size(); p++){
                    if(protDB.get(p).get(0).equals(locus))
                        orthologDB.add(start, protDB.get(p));
                }
                matchesNum++;

                //add to locusList
                end = start+matchesNum-1;
                locusList.add(new ArrayList(Arrays.asList(locus, Integer.toString(matchesNum), Integer.toString(start), Integer.toString(end))));
                start = end+1;

                Utilities.setTaskProgress(++i);
            }
        }
    }

    /**
     * creating alignment database
     * using ClustalW2
     */
    private void createSeqAlignment() throws Exception {

        Utilities.startTask("Task 3 of 4: Performing Sequence Alignment", orthologDB.size());

        if(locusSList.isEmpty())
            return;

        String locus = "";
        String HID = "";
        boolean newSet = true;
        BufferedWriter inputFile = null;
        int locusIndex = 0;
        int start = Integer.valueOf(locusList.get(locusIndex).get(2));
        int end = Integer.valueOf(locusList.get(locusIndex).get(3));

        //locusList contains the start (index=2) and end (index=3) indices of orthologDB for each ortholog group
        for (int i = 0; i < orthologDB.size(); i++){
            //only do alignment if there is at least two proteins in this ortholog group
            if(start != end){
                if(i == start){
                    locus = orthologDB.get(i).get(0);
                    HID = Utilities.getHID(locus);
                    //create input file
                    inputFile = new BufferedWriter(new FileWriter(alnOutFolderPath+HID+".fasta"));
                }

                inputFile.write(">" + Utilities.getHIDSpecies(orthologDB.get(i).get(1)) + "\n" + orthologDB.get(i).get(4) + "\n");

                if(i == end){
                    if(inputFile != null)
                        inputFile.close();
                    
                    /* 
                    //run local ClustalW2 executable
                    String ClustalCommand = "\""+exeFolderPath+"clustalw2.exe\" -infile="+"\""+alnOutFolderPath+HID+".fasta\" -outfile="+"\""+alnOutFolderPath+HID+".aln\"";
                    ClustalCommand += " -align -type=protein -outorder=input -seqnos=on";
                    if(!Utilities.sendToExec(ClustalCommand, false))
                        Utilities.displayMessage("An error occurred while trying to execute ClustalW2 for "+locus, parentFrame);
                    ClustalCommand = "\""+exeFolderPath+"clustalw2.exe\" -infile="+"\""+alnOutFolderPath+HID+".aln\" -tree";
                    if(!Utilities.sendToExec(ClustalCommand, false))
                        Utilities.displayMessage("An error occurred while trying to execute ClustalW2 for "+locus, parentFrame);
                    */
                    
                    //run local Muscle executable
                    //commented out exeFolderPath, to make it work on Linux - muscle is included in path variables
                    //to make it work on windows again, add muscle to path or uncomment the next line (and comment the linux-specific command
                    
                    //calling muscle on a path
                    String[] ClustalCommand = {exeFolderPath+"muscle", "-in", alnOutFolderPath+HID+".fasta", "-out", alnOutFolderPath+HID+".aln", "-tree1", alnOutFolderPath+HID+".ph", "-maxiters", "2", "-clwstrict"};
                    
                    //calling muscle, assuming program is defined in PATH
                    //String ClustalCommand = "muscle -in "+""+alnOutFolderPath+HID+".fasta -out "+""+alnOutFolderPath+HID+".aln -tree1 "+alnOutFolderPath+HID+".ph -maxiters 2 -clwstrict";
                    
                    if(!Utilities.sendToExec(ClustalCommand, false))
                        Utilities.displayMessage("An error occurred while trying to execute Muscle for "+locus, parentFrame);
                }
            }
            
            if(i == end){
                locusIndex++;
                if(locusIndex < locusList.size()){
                    start = Integer.valueOf(locusList.get(locusIndex).get(2));
                    end = Integer.valueOf(locusList.get(locusIndex).get(3));
                }
            }

            Utilities.setTaskProgress(i+1);
        }
    }

    /**
     * generate conservation score based on algorithm derived from CI value calculations
     * @modifies cscoresDB
     * @requires output from the .aln alignment file must list the input species first, before all other species
     */
    private void createCScoresDB() throws Exception {
        createPhosphositeDB();
        Utilities.startTask("Task 4 of 4: Calculating Conservation Scores", (int) (phosphositeDB.size()+phosphositeDB.size()*0.1));
        
        if(locusSList.isEmpty()){
            Utilities.setTaskMax(100);
            Utilities.setTaskProgress(100);
            return;
        }

        //parse .aln file into the following data structure:
        //2D array seqAlignment: [[aa, aa, aa, etc.]]
        //with index of first level representing residue pos (with respect to original species)
        //and second level representing all the aligned residues at that pos
        ArrayList<ArrayList<Character>> seqAlignment = new ArrayList<ArrayList<Character>>();
        ArrayList<ArrayList<Integer>> gapsN = new  ArrayList<ArrayList<Integer>>();
        short orthologsN = 1; //keeps track of number of orthologs for each locus, easier to count than to find it in locusList

        ArrayList<String> speciesList = null;
        int i = 0;
        String prevLocus = "";
        
        for(ArrayList<String> phosphosite : phosphositeDB){
            
            String locus = phosphosite.get(0);
            String HID = Utilities.getHID(locus);
            String locus_species = phosphosite.get(1);
            String inputSpecies = Utilities.getSpecies(locus_species, false);
            char siteAA = phosphosite.get(2).charAt(0);

            //phosphosite location (for input species only); pos relative to protein seq (count starts at 1)
            int siteN = Integer.parseInt(phosphosite.get(2).substring(1));
            
            if (!locus.equals(prevLocus)) {
                prevLocus = locus;
                speciesList = new ArrayList<String>();
                
                //sequence alignment
                //each row (first level of 2D array) represents residue pos
                //each column (second level of 2D array) represents the different orthologs
                //pos starts with 0
                seqAlignment = new ArrayList<ArrayList<Character>>();

                //counts for the gaps
                //each row represents the different orthologs; each column represents the pos where the gap presents
                //pos starts with 0
                gapsN = new  ArrayList<ArrayList<Integer>>();

                int lastMax = 0; //tracks the starting index for other species for each block of alignments
                short setNewN = 0; //setNew occurrences
                short inputSpeciesIdx = -1;
                boolean firstTime = true;
                String firstSpecies = null; //first species in the order of species output .aln list
                
                BufferedReader inputStream = null;
                String s;
                try {
                    inputStream = new BufferedReader(new FileReader(alnOutFolderPath + HID + ".aln"));
                    while ((s = inputStream.readLine()) != null) {
                        if (s.startsWith(HID)) {
                            String sCleaned = s.substring(s.indexOf(' '));
                            char[] seq = sCleaned.toCharArray();
                            
                            String species = s.substring(0,s.indexOf(' '));
                            species = species.replace('_', ' ').substring(species.indexOf('|')+1);
                            //NOTE: some alignment software truncates the name in the output, so 
                            //species name comparison will be based startsWith instead of equals
                            
                            if(firstTime == true){
                                firstSpecies = species;
                                firstTime = false;
                            }
                            
                            boolean isFirstSpecies = firstSpecies.startsWith(species) ? true : false;
                            boolean isInputSpecies = inputSpecies.startsWith(species) ? true : false;
                            
                            if (isFirstSpecies) {
                                lastMax = seqAlignment.size();
                                orthologsN = 1;
                                setNewN++;
                            } else
                                orthologsN++;

                            if(inputSpecies.startsWith(species))
                                inputSpeciesIdx = (short) (orthologsN - 1);
                            
                            int pos = isFirstSpecies ? 0 : lastMax; //the 0 here is random; pos is useful when setNew is false

                            if(setNewN == 1){
                                gapsN.add(new ArrayList<Integer>());
                                if(isInputSpecies)
                                    speciesList.add(0, species);
                                else
                                    speciesList.add(species);
                            }

                            int a = 0;
                            for (int j = 0; j < seq.length; j++) {
                                if (Character.isLetter(seq[j]) || seq[j] == '-') {
                                    if (isFirstSpecies) {
                                        seqAlignment.add(new ArrayList<Character>(Arrays.asList(seq[j])));
                                    } else if (isInputSpecies){
                                        seqAlignment.get(pos++).add(0, seq[j]);
                                    } else {
                                        seqAlignment.get(pos++).add(seq[j]);
                                    }

                                    if(seq[j] == '-')
                                        gapsN.get(orthologsN-1).add(lastMax+a); 
                                    //Note: gapsN is not the same order as seqAlignment and speciesList
                                    //will fix this later. Easier to save entries in gapsN in the same
                                    //order as that listed in alignment file for now.

                                    a++;
                                }
                            }
                        }
                    }

                    if (inputStream != null) {
                        inputStream.close();
                    }
                } catch (FileNotFoundException e) {
                    continue;
                }
                
                
                //move the entries for input species to index 0
                if(inputSpeciesIdx != 0){
                    gapsN.add(0,gapsN.get(inputSpeciesIdx));
                    gapsN.remove(inputSpeciesIdx+1);
                }
                //seqAlignment and speciesList are already in this order
                        
            }

            //tracks the phosphosite location of each species
            //pos relative to protein sequence
            int[] sitesN = new int[orthologsN];
            sitesN[0] = siteN;
            
            //get site location with ('-') in the alignment taken into account (for input species)
            //siteNabs starts at 1 but takes gaps into account
            int siteNabs = siteN;
            for(int y = 0; y < gapsN.get(0).size();y++){
                if(gapsN.get(0).get(y) < siteN+y)
                    siteNabs++;
            }

            //for the other species, start with the site = siteNabs
            for(int s = 1; s < orthologsN; s++)
                sitesN[s] = siteNabs;

            //find the phosphosite location of other species
            for(int x = 1; x < gapsN.size(); x++){
                for(int y = 0; y < gapsN.get(x).size();y++){
                    if(gapsN.get(x).get(y) < siteNabs)
                        sitesN[x]--;
                }
            }

            String phosphositesList = "";
            //put list of corresponding phosphosite in other species together and add to phosphositeDB
            //format: species:phosphosite;
            for(int s = 0; s < orthologsN; s++){
                phosphositesList += speciesList.get(s)+":";
                String letter = String.valueOf(seqAlignment.get(siteNabs-1).get(s));
                phosphositesList += (letter.equals("-")) ? letter : letter+sitesN[s];
                phosphositesList += "; ";
            }
            phosphositeDB.get(i).add(phosphositesList);

            //now get the specified width of sequence residues surrouding each phosphosite
            //if phosphosite is positioned too close to C-termious or N-termious such that
            //the width of site alignment cannot be complete, then add blanks until width is met
            //currently the blanks are not taken into account for conservation scoring
            //additional feature can be implemented such as the blanks can be replaced with letters
            
            int startOffset = 0, endOffset = 0;
            int start = siteNabs-siteAlignNgbrLen-1;
            int end = siteNabs+siteAlignNgbrLen-1;
            int endOfSeqAlign = seqAlignment.size()-1;
            if(start < 0)
                startOffset = -1*start;
            if(end > endOfSeqAlign)
                endOffset = end - endOfSeqAlign;

            //grab the specified width of sequence around phophosite
            phosphositeAlignDB.add(new char[siteAlignTotalLen][orthologsN]);
            for(int z = 0; z < siteAlignTotalLen; z++){
                if(z < startOffset){
                    for(int g = 0; g < orthologsN; g++)
                        phosphositeAlignDB.get(i)[z][g] = ' ';
                } else if(z >= siteAlignTotalLen-endOffset){
                    for(int g = 0; g < orthologsN; g++)
                        phosphositeAlignDB.get(i)[z][g] = ' ';
                } else {
                    for(int g = 0; g < orthologsN; g++)
                        phosphositeAlignDB.get(i)[z][g] = seqAlignment.get(start+z).get(g);
                }
            }

            Utilities.setTaskProgress(++i);
        }

        //calculate conservation score
        //this class will take the phosphositeAlignDB and fill in the scores to cscoresDB
        ConservationScore cs = new ConservationScore(this, false);
        Utilities.setTaskProgress((int) (phosphositeDB.size()+phosphositeDB.size()*0.1));
    }
}
