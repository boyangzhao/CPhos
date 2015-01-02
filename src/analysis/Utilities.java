package analysis;

import java.io.*;
import java.net.*;
import java.util.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import gui.*;

/**
 * @author Boyang Zhao
 * @description set of utility methods to be used by analysis and gui classes
 */
public class Utilities {
    
    /* basic operations*/
    
    public static double log(double val, double base){
        return (Math.log(val)/Math.log(base));
    }
    
    public static String reportRunningTime(long runningtime){
        float seconds = (float)runningtime/(float)1000.0;
        
        if(seconds < 60)
            return String.format("%.2f",seconds)+" seconds";
        else if(seconds < 3600)
            return String.format("%.2f",seconds/60)+" minutes";
        else
            return String.format("%.2f",seconds/3600)+" hours";
    }
    
    public static String getRuntimeStats(long runningtime){
        String msg =  "Total running time: " + runningtime + " millis\n";
        msg  += "Total running time: " + Utilities.reportRunningTime(runningtime) + "\n";
        msg  += "Total available memory: " + Runtime.getRuntime().totalMemory()/(1024*1024) + " MB\n";
        msg  += "Used memory: " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/(1024*1024) + " MB\n";
        msg  += "Free memory: " + Runtime.getRuntime().freeMemory()/(1024*1024) + " MB\n";
        msg  += "Maximum available memory: " + Runtime.getRuntime().maxMemory()/(1024*1024) + " MB";
        
        return msg;
    }

    /**
     * extract gene symbols or protein accession number list from hashset
     * the hashset is from object in serialized homologene database
     * @param set hashset containing list of genes (each entry is one unique gene) and one entry of a proteins ("||" as delimiter)
     * @param genes look below (@return)
     * @return returns gene symbols list if true, otherwise will return protein accession number list
     */
    public static String getFromSet(HashSet<String> set, boolean genes){
        String joinedgenes = "";
        Iterator<String> iter = set.iterator();
        while(iter.hasNext()){
            String tmp = iter.next();
            if(tmp.contains("||")){
                if(!genes) //can just return now, don't need to contine to get all the genes
                    return tmp;
            } else
                joinedgenes += tmp+",";
        }
        
        if(!genes)
            return null;
        else
            return joinedgenes.substring(0,joinedgenes.lastIndexOf(','));
    }

    /* updating methods*/
    public static void startTask(String taskname) {
        if(Analysis.display){
            CPhos.taskProgressBar.setIndeterminate(true);
            CPhos.progressTaskLabel.setText(taskname);
        }
    }

    public static void startTask(String taskname, int max) {
        if(Analysis.display){
            CPhos.taskProgressBar.setMaximum(max);
            CPhos.taskProgressBar.setValue(0);
            CPhos.taskProgressBar.setIndeterminate(false);
            CPhos.progressTaskLabel.setText(taskname);
        }
    }

    public static void setOverallMax(int val) {
        if(Analysis.display)
            CPhos.overallProgressBar.setMaximum(val);
    }

    public static void setOverallProgress(int val) {
        if(Analysis.display){
            CPhos.overallProgressBar.setValue(val);
            switch (val/25) {
                case 1:
                    CPhos.overviewLabel1.setFont(CPhos.defaultfontBold);
                    break;
                case 2:
                    CPhos.overviewLabel2.setFont(CPhos.defaultfontBold);
                    break;
                case 3:
                    CPhos.overviewLabel3.setFont(CPhos.defaultfontBold);
                    break;
                case 4:
                    CPhos.overviewLabel4.setFont(CPhos.defaultfontBold);
                    break;
            }
        }
    }

    public static void setTaskMax(int val) {
        if(Analysis.display)
            CPhos.taskProgressBar.setMaximum(val);
    }
    
    public static void setTaskProgress(int val) {
        if(Analysis.display)
            CPhos.taskProgressBar.setValue(val);
    }

    public static void displayMessage(String msg, java.awt.Window parent) {
        if(Analysis.display)
            JOptionPane.showMessageDialog(parent, msg);
        else
            System.out.println(msg);
    }

    /**
     * update settings: update all values of settings variables (all static)
     * in Analysis and ConservationScore class
     * @modifies static variables in Analysis and ConservationScore class
     */
    public static void updateSettings(){
        try{
            String userdir = System.getProperty("user.dir")+System.getProperty("file.separator");
            String slash = System.getProperty("file.separator");
            
            Properties p = new Properties();
            p.load(new FileInputStream("./settings.ini"));
            Analysis.dbFolderPath = p.getProperty("dir.database");
            Analysis.exeFolderPath = p.getProperty("dir.executable");
            Analysis.alnOutFolderPath = p.getProperty("dir.output");
            if(Analysis.dbFolderPath.isEmpty()) Analysis.dbFolderPath = userdir+"db"+slash;
            if(Analysis.exeFolderPath.isEmpty()) Analysis.exeFolderPath = userdir+"exe"+slash;
            if(Analysis.alnOutFolderPath.isEmpty()) Analysis.alnOutFolderPath = userdir+"outputs"+slash;
            
            
            Analysis.db_species = new ArrayList(Arrays.asList(p.getProperty("database.species").split("\\|\\|")));
            Analysis.db_target = new ArrayList(Arrays.asList(p.getProperty("database.target").split("\\|\\|")));
            Analysis.siteAlignNgbrLen = Integer.parseInt(p.getProperty("as.siteAlignmentNgbrLenOneSide"));
            Analysis.cPlaceholder = p.getProperty("as.cPlaceHolder").charAt(0);
            Analysis.nPlaceholder = p.getProperty("as.nPlaceHolder").charAt(0);
            ConservationScore.logbase = Double.parseDouble(p.getProperty("cs.logarithmBase"));

            //TODO: finish default
            //ConservationScore.bgfreq_a = new ArrayList(Arrays.asList(0.05,0.05,0.01,0.09,0.05,0.05,0.05,0.05,0.01,0.01,0.13,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
            ConservationScore.bgfreq_a = new ArrayList(Arrays.asList(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
        } catch (Exception e) {
            displayMessage("Unable to find settings.ini! Default settings will be used.", null);
            
            String userdir = System.getProperty("user.dir")+System.getProperty("file.separator");
            String slash = System.getProperty("file.separator");
            
            Analysis.dbFolderPath = userdir+"db"+slash;
            Analysis.exeFolderPath = userdir+"exe"+slash;
            Analysis.alnOutFolderPath = userdir+"outputs"+slash;
            Analysis.siteAlignNgbrLen = 6;
            Analysis.cPlaceholder = 'B';
            Analysis.nPlaceholder = 'Z';
            ConservationScore.logbase = 2.0;

            //TODO: what should be the default
            ConservationScore.bgfreq_a = new ArrayList(Arrays.asList(0.05,0.05,0.01,0.09,0.05,0.05,0.05,0.05,0.01,0.01,0.13,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));

            //list all .sdb as possible databases
            Analysis.db_species = new ArrayList<String>();
            Analysis.db_target = new ArrayList<String>();
            File dbDir = new File(Analysis.dbFolderPath);
            String[] dbs = dbDir.list();
            if (dbs != null) {
                for (int i = 0; i < dbs.length; i++) {
                    if (dbs[i].endsWith(".sdb") && !dbs[i].startsWith("homologene.")) {
                        Analysis.db_species.add(dbs[i]);
                        Analysis.db_target.add(dbs[i]);
                    }
                }
            }
        }
    }
    
    /**
     * write to settings.ini
     * @param propertyname propertyname to update
     * @param propertyval property add to append, delimited by ||
     */
    public static void writeToSettings(String propertyname, String propertyval){
        Properties p = null;
        try {
            p = new Properties();
            p.load(new FileInputStream("./settings.ini"));
            String dbspecies = p.getProperty(propertyname);
            p.setProperty(propertyname, propertyval);
            p.store(new BufferedOutputStream(new FileOutputStream("settings.ini")), "default settings for CPhos\n\nauthor: Boyang Zhao" +
                    "\n\ninstructions: CPhos looks for settings.ini in the same directory as program (.jar) to initialize setting variables" +
                    "\nif this file doesn't exist, the program will start with the default values\n\nLast modified: ");
        } catch(Exception e){}
    }

    /* file handling methods */
    /**
     * get FASTA with the input protein id from NCBI Web Service (EUtilities)
     * @param proteinID protein refseq accession number
     * @param taxonomyclass target database (aka the taxonomy class to be limited to)
     * @return string array [found, species, accession number, gene symbol, protein seq]
     *         value of index 1 in return array is "found" is the protein is classified in the input taxonomy class; "notfound" otherwise
     */
    public static String[] callws_getgp(String proteinID, String taxonomyclass){
        String[] parsed = new String[5];
        parsed[0] = "notfound";
        boolean success = false;
        while (!success) {
            try {
                URL ncbi = new URL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=" + proteinID + "&rettype=gp&retmode=text"); //Call NCBI protein
                BufferedReader in = new BufferedReader(
                                    new InputStreamReader(
                                    ncbi.openStream()));
                String s = null;
                while ((s = in.readLine()) != null) {
                    if (s.startsWith("  ORGANISM")) {
                        //species
                        parsed[1] = s.replace("  ORGANISM  ", "");

                        //continue search to see if the input protein is in the taxonomy class
                        while((s = in.readLine()) != null){
                            //if out of the organism section, then break loop
                            if(!s.startsWith("   "))
                                break;
                            else {
                                String search = s.toLowerCase();
                                if(search.contains(taxonomyclass))
                                    parsed[0] = "found";
                            }
                        }
                    } else if (s.startsWith("LOCUS")){
                        //accession number
                        String[] tmp = s.split("\\s+");
                        parsed[2] = tmp[1];
                    } else if (s.contains("/gene=")){
                        //gene symbol
                        parsed[3] = s.trim().replace("/gene=", "").replace("\"", "");
                    } else if (s.startsWith("ORIGIN")){
                        //protein sequence
                        String sequence = "";
                        do {
                            if(s.contains("//"))
                                break;
                            
                            if (!s.startsWith("ORIGIN")) {
                                sequence += s.replaceAll("[\\s\\d]+", "");
                            }
                        } while ((s = in.readLine()) != null);
                        parsed[4] = sequence.toString().toUpperCase();
                    }
                }
                in.close();
                success = true;
            } catch (MalformedURLException e) {
            } catch (IOException e) {
                success = false;
            }
        }
        return parsed;
    }
    
    /**
     * deserialize 2D array object (refseq databases)
     * @param filepath filepath of database
     * @param buffersize buffer size
     * @return 2D array of the deserialized object
     */
    public static ArrayList<ArrayList<String>> deserializeDB(String filepath, int buffersize) throws Exception {
        ObjectInputStream inputStream = null;
        inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(filepath), buffersize));
        ArrayList<ArrayList<String>> tmp = (ArrayList<ArrayList<String>>) inputStream.readObject();
        if (inputStream != null){
            inputStream.close();
            return tmp;
        }
        
        return null;
    }
    
    /**
     * create the .sdb protein database using .fasta of target database and homologene data
     * this .sdb file consists of one object
     *  1) ArrayList<ArrayList<String>> that is 2D array of the data
     *          format: [[locus_species, protein accession number, gene symbol, protein sequence]]
     *          note: locus_species has the format: HID|gene symbols[species]
     * @param dbFilePath file path of .fasta file
     * @param dbname database name
     * @return true if successful, false otherwise
     */
    public static boolean serializeDB(String dbFilePath, String dbname, String destFoldername, javax.swing.JProgressBar dbProgressBar){
        //get the most recent version of homologene database (.sdb version)
        String hfilename = gethDBfilename();

        ObjectInputStream inputStream = null;
        ArrayList<ArrayList<String>> hDB = new ArrayList<ArrayList<String>>();
        HashMap<String, HashSet> hDB_hash = new HashMap<String, HashSet>();
        try {
            inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(Analysis.dbFolderPath+hfilename), 200000));
            hDB = (ArrayList<ArrayList<String>>) inputStream.readObject();
            hDB_hash = (HashMap<String, HashSet>) inputStream.readObject();
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        } catch (ClassNotFoundException e) {
        } finally {
            try {
                if (inputStream != null)
                    inputStream.close();
            } catch (IOException e) {}
        }

        //set progress bar
        File f = new File(dbFilePath);
        long size = f.length();
        dbProgressBar.setMaximum(100);
        
        //parse fasta file into hash table
        BufferedReader fastaOutputStream = null;
        //key = refseqID (without version number); value = [species]protein sequence
        HashMap<String, String> protseq = new HashMap<String, String>();
        //key = refseqID (wihtout version number); value = refseqID (with version number)
        HashMap<String, String> protID = new HashMap<String, String>();
        try {
            fastaOutputStream = new BufferedReader(new FileReader(dbFilePath));
            FileReader a = new FileReader(dbFilePath);
            
            //parse fasta file first into array
            //fasta first line has the format: >gi|gi number|ref|refseq number|protein name [species]
            //note: it's possible that the protein name also contains '[' and ']'
            String seq = ""; //value format: [species]protein sequence
            String ref = ""; //key format: refseq number
            String s;
            long csize = 0;
            Boolean skip = false;
            while ((s = fastaOutputStream.readLine()) != null) {
                if(s.startsWith(">")){
                    //hash the previous seq and start over
                    if (!seq.isEmpty()){
                        seq = seq.trim();
                        //remove the version number
                        String refwov = ref.substring(0, ref.lastIndexOf('.'));
                        protseq.put(refwov, seq);
                        protID.put(refwov, ref);
                    }
                    
                    //start collecting new seq
                    int lindex = s.lastIndexOf('[');
                    if(lindex == -1){
                        Utilities.displayMessage("The protein name in \""+s+"\" does not contain [species]. This entry will be skipped", null);
                        skip = true;
                    } else {
                        seq = s.substring(lindex);
                        ref = s.substring(s.indexOf("ref|")+4, s.indexOf('|',s.indexOf("ref|")+4));
                        skip = false;
                    }
                } else if(s.length() > 0) {
                    if(!skip){
                        if(Character.isLetter(s.charAt(0))){
                            seq += s;
                        }
                    } else
                        seq = "";
                }

                csize += s.getBytes().length;
                dbProgressBar.setValue((int)(csize*100/(2*size)));
            }

            //hash the last one
            seq = seq.trim();
            String refwov = ref.substring(0, ref.lastIndexOf('.'));
            protseq.put(refwov, seq);
            protID.put(refwov, ref);
        } catch (FileNotFoundException e) {
            return false;
        } catch (IOException e) {
            return false;
        } finally {
            try {
                if (fastaOutputStream != null) {
                    fastaOutputStream.close();
                }
            } catch (IOException e) {
                return false;
            }
        }

        //create .sdb database
        ArrayList<ArrayList<String>> sDB = new ArrayList<ArrayList<String>>();
        String sdbFileName = createDBFileName(dbFilePath, dbname);
        if(destFoldername.endsWith("\\") || destFoldername.endsWith("/"))
            destFoldername += System.getProperty("file.separator");
        BufferedWriter txtOutputStream = null;
        ObjectOutputStream objOutputStream = null;
        try{
            //write the .sdb file
            txtOutputStream = new BufferedWriter(new FileWriter(destFoldername+sdbFileName+".txt", false));
            objOutputStream = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(destFoldername+sdbFileName, false)));

            Iterator<ArrayList<String>> iter = hDB.iterator();
            int i = 0;
            while(iter.hasNext()){
                ArrayList<String> item = iter.next();
                //if the current protein ID in homologene database (item.get(5)) is also in the input species fasta database (proseq key)
                //then continue to put this in the .sdb database
                //the value of proseq has the format: [species]protein sequence
                if(protseq.containsKey(item.get(5))){
                    String speciesseq = protseq.get(item.get(5));
                    String species = speciesseq.substring(0,speciesseq.indexOf(']')+1);
                    String seq = speciesseq.substring(speciesseq.indexOf(']')+1);
                    String[] tmp = new String[4];
                    String genesymbols = getFromSet(hDB_hash.get(item.get(0)), true);
                    tmp[0] = item.get(0)+"|"+genesymbols+"_"+species;
                    tmp[1] = protID.get(item.get(5));
                    tmp[2] = item.get(3);
                    tmp[3] = seq;
                    sDB.add(new ArrayList<String>(Arrays.asList(tmp)));
                    txtOutputStream.append("["+tmp[0]+","+tmp[1]+","+tmp[2]+","+tmp[3]+"]\n");
                }

                dbProgressBar.setValue((int)((0.5+(0.5)*(++i)/hDB.size())*100));
            }
            
            objOutputStream.writeObject(sDB);
            
            //update/write settings.ini file
            Properties p = new Properties();
            p.load(new FileInputStream("./settings.ini"));
            String dbspecies = p.getProperty("database.species");
            writeToSettings("database.species",dbspecies+"||"+sdbFileName);
        } catch (FileNotFoundException e) {
            return false;
        } catch (IOException e) {
            return false;
        } finally {
            try {
                if (txtOutputStream != null) {
                    txtOutputStream.close();
                }
                if (objOutputStream != null) {
                    objOutputStream.close();
                }
            } catch (IOException e) {
                return false;
            }
        }
        
        dbProgressBar.setStringPainted(true);
        dbProgressBar.setString("Complete");
        
        return true;
    }

    /**
     * serialize the .data file and create the .sdb file (homologene database)
     * entries where there are multiple proteins from the same species and same HID are removed
     * .data file format (tab delimited): HID, Taxonomy ID, Gene ID, Gene Symbol, Protein GI, Protein Accession Number
     * this .sdb file consists of two objects
     *  1) ArrayList<ArrayList<String>> that is 2D array of the data
     *          format: [[HID, Taxonomy ID, Gene ID, Gene Symbol, Protein GI, Protein Accession Number]]
     *  2) HashMap<String, HashSet> with key=HID and value=hashset of gene symbols (one gene symbol per entry)
     *                                                     and one additional entry of all proteins in this HID group ("||" as delimiter)
     * @param dbFilePath file path of .data file
     * @param dbversion build version of databae
     * @return true if successful, false otherwise
     */
    public static boolean serializeDB(String dbFilePath, Integer dbbuild, javax.swing.JProgressBar dbProgressBar){
        //set progress bar
        File f = new File(dbFilePath);
        long size = f.length();
        dbProgressBar.setMaximum(100);

        String sdbFileName = Analysis.dbFolderPath+"homologene."+dbbuild.toString()+".sdb";
        ArrayList<ArrayList<String>> hDB = new ArrayList<ArrayList<String>>();
        //for hDB_hash: key=HID, value=gene symbols (+1 entry with list of proteins)
        HashMap<String, HashSet> hDB_hash = new HashMap<String, HashSet>();
        BufferedReader inputStream = null;
        BufferedWriter txtOutputStream = null;
        ObjectOutputStream objOutputStream = null;
        try {
            long csize = 0;
            String s;
            TreeMap<String, String> inputData = new TreeMap<String, String>();
            int max = 0;
            inputStream = new BufferedReader(new FileReader(dbFilePath));
            txtOutputStream = new BufferedWriter(new FileWriter(sdbFileName+".txt", false));
            objOutputStream = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(sdbFileName, false)));

            //read in the homologene data
            //if there are more than one proteins that share same HID and taxonomy ID, max the value with * instead
            //later ones with * value will be omitted
            while ((s = inputStream.readLine()) != null) {
                ArrayList<String> tmp = new ArrayList(Arrays.asList(s.split("\\t")));
                //key here is HID\ttaxonomyID
                String key = tmp.get(0)+"\t"+tmp.get(1);
                Boolean exist = inputData.containsKey(key);
                inputData.put(key, exist ? "*" : s);

                csize += s.getBytes().length;
                dbProgressBar.setValue((int)(csize*100/(2*size)));
            }

            Iterator iter = inputData.entrySet().iterator();
            int c = 0;
            while(iter.hasNext()){
                Map.Entry entry = (Map.Entry) iter.next();
                c++;

                if(entry.getValue() != "*"){
                    ArrayList<String> row = new ArrayList(Arrays.asList(entry.getValue().toString().split("\\t")));

                    //homologene has protein that contains version number
                    //get rid of this first; we'll get the most recent version of protein info later
                    row.set(5, row.get(5).substring(0, row.get(5).lastIndexOf('.')));
                    
                    //remove any '[' and ']' in the gene name
                    row.set(3, row.get(3).replace("[",""));
                    row.set(3, row.get(3).replace("]",""));

                    hDB.add(row);
                    
                    if(hDB_hash.containsKey(row.get(0))){
                        //this entry has the same HID as previous
                        hDB_hash.get(row.get(0)).add(row.get(3).toLowerCase());

                        String protlist = getFromSet(hDB_hash.get(row.get(0)),false);
                        hDB_hash.get(row.get(0)).remove(protlist);
                        hDB_hash.get(row.get(0)).add(protlist+row.get(5)+"||");
                    } else {
                        //new HID
                        HashSet<String> gset = new HashSet<String>();
                        gset.add(row.get(3).toLowerCase());
                        gset.add(row.get(5)+"||");
                        hDB_hash.put(row.get(0), gset);
                    }
                }

                dbProgressBar.setValue((int)(50+(c*100/(2*inputData.size()))));
            }
            
            //remove any entries where there is only one protein in the HID group
            Iterator iter2 = hDB.iterator();
            while(iter2.hasNext()){
                ArrayList<String> entry = (ArrayList<String>) iter2.next();
                String HIDkey = entry.get(0);
                
                String protlist = getFromSet(hDB_hash.get(HIDkey),false);
                String[] proteins = protlist.split("\\|+");
                if(proteins.length < 2){
                    iter2.remove();
                } else {
                    txtOutputStream.append("["+entry.get(0)+","+entry.get(1)+","+entry.get(2)+","+entry.get(3)+","+entry.get(4)+","+entry.get(5)+"]\n");
                }
            }
            
            objOutputStream.writeObject(hDB);
            objOutputStream.writeObject(hDB_hash);
        } catch (FileNotFoundException e) {
            return false;
        } catch (IOException e) {
            return false;
        } finally {
            try {
                if (inputStream != null) {
                    inputStream.close();
                }
                if (txtOutputStream != null) {
                    txtOutputStream.close();
                }
                if (objOutputStream != null) {
                    objOutputStream.close();
                }
            } catch (IOException e) {
                return false;
            }
        }

        dbProgressBar.setStringPainted(true);
        dbProgressBar.setString("Complete");
        
        return true;
    }
    
    /**
     * JFileChooser for saving file, send it to main saveFile method
     */
    public static String saveFile(java.awt.Window parent, String filename) {
        try{
            return saveFilePathGiven(parent, new File(filename).getCanonicalPath());
        } catch(IOException e){}

        return "";
    }

    /**
     * JFileChooser for saving file; has the ability to check if the file already exists
     * @param parent parent window
     * @param path initial path for where the file is to be saved to
     * @return the final chosen path for where the file is to be saved to
     */
    public static String saveFilePathGiven(java.awt.Window parent, String path) {
        try {
            JFileChooser dbdirSelect = new JFileChooser();
            dbdirSelect.setSelectedFile(new File(path));
            String dbFilePath = "";
            int returnVal = dbdirSelect.showDialog(parent, "Save File");
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File dbFile = dbdirSelect.getSelectedFile();
                dbFilePath = dbFile.getPath();

                if (dbFile.exists()) {
                    int val = JOptionPane.showConfirmDialog(parent,
                            "File already exists! Do you wish to replace existing file?",
                            "Warning",
                            JOptionPane.YES_NO_OPTION);
                    if(val == JOptionPane.YES_OPTION)
                        return dbFilePath;
                    else
                        return saveFile(parent, dbdirSelect.getSelectedFile().getCanonicalPath());
                }

                return dbFilePath;
            }
        } catch(IOException e){}

        return "";
    }

    /**
     * JFileChooser for selecting file
     * @param extRestrict
     * @param errorMsg
     * @param parent
     * @return
     */
    public static String selectFile(String extRestrict, String errorMsg, java.awt.Window parent) {
        JFileChooser dbdirSelect = new JFileChooser();
        String dbFilePath = "";
        int returnVal = dbdirSelect.showDialog(parent, "Select File");
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File dbFile = dbdirSelect.getSelectedFile();
            dbFilePath = dbFile.getPath();

            if (!dbFilePath.endsWith(extRestrict)) {
                Utilities.displayMessage(errorMsg, parent);
                return "";
            } else {
                return dbFilePath;
            }
        }

        return "";
    }

    /**
     * JFileChooser for selecting folder
     * @param parent
     * @return
     */
    public static String selectFolder(java.awt.Window parent) {
        JFileChooser dbdirSelect = new JFileChooser();
        dbdirSelect.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        dbdirSelect.setAcceptAllFileFilterUsed(false);
        int returnVal = dbdirSelect.showDialog(parent, "Select Folder");

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File fldr = dbdirSelect.getSelectedFile();
            String dbFolderPath = fldr.getPath();

            return dbFolderPath;
        }

        return "";
    }

     /**
     * open documentation file with system default reader
     * @param filename filepath
     */
    public static void openDocument(final String filepath, final java.awt.Window window){
        new Thread() {
            @Override
            public void run() {
                try {
                    URI uri = URI.create(filepath);
                    java.awt.Desktop desktop;
                    if (java.awt.Desktop.isDesktopSupported()) {
                        desktop = java.awt.Desktop.getDesktop();
                        desktop.browse(uri);
                    } else {
                        String errorMsg = "Unable to open documentation. Please browse to " + filepath + " to open it manually.";
                        Utilities.displayMessage(errorMsg, window);
                    }
                } catch (IOException e) {
                    Utilities.displayMessage("Unable to open/find documentation. Please try to look for " + filepath + " manually.", window);
                }
            }
        }.start();
    }
    
    public static String readFile(String filename){
        String resultString = "";
        BufferedReader f = null;
        try {
            f = new BufferedReader(new FileReader(filename));
            String s;
            boolean beg = true;
            while ((s = f.readLine()) != null) {
                if(beg){
                    resultString = s;
                    beg = false;
                } else
                    resultString = resultString+"\n"+s;
            }
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        } finally {
            try {
                if (f != null)
                    f.close();
            } catch(IOException e){}
        }
        
        return resultString;
    }
    
    /**
     * creates the .sdb filename using input dbname and the date of modification of fasta database (located in dbfilepath)
     * e.g. input: Rattus norvegicus
     * e.g. output: rattus_norvegicus_[refseq].07_01_10.sdb
     */
    public static String createDBFileName(String dbfilepath, String dbname){
        File f = new File(dbfilepath);
        long lastmod = f.lastModified();
        DateFormat formatter = new SimpleDateFormat("MM_dd_yy");
        Calendar cal = Calendar.getInstance();
        cal.setTimeInMillis(lastmod);
        String dbfilename = dbname.toLowerCase().replace(' ', '_');
        dbfilename += "_[refseq]."+formatter.format(cal.getTime())+".sdb";
        return dbfilename;
    }

    /**
     * creates the web service filename using input dbname
     * e.g. input: Mammalia
     * e.g. output: mammalia_[refseq web service].ws
     */
    public static String createDBFileName(String dbname){
        return dbname.toLowerCase().replace(' ', '_')+"_[refseq_web_service].ws";
    }
    
    /**
     * extract the database name (format: Species [refseq]) from the database file name (format: species_[refseq].mmddyyy.sdb)
     * -or-
     * extract the database name (format: Species [refseq web service]) from the database file name (format: species_[refseq_web_service].ws)
     *
     * e.g. input: rattus_norvegicus_[refseq].07_01_10.sdb
     * e.g. output if withdbsource is true: Rattus norvegicus [refseq]
     * e.g. output if withdbsource is false: Rattus norvegicus
     */
    public static String getDBName(String dbfilename, Boolean withdbsource){
        String tmp = "";
        if(dbfilename.endsWith(".sdb"))
            tmp = dbfilename.substring(0, dbfilename.indexOf(".sdb")-9);
        else if(dbfilename.endsWith(".ws"))
            tmp = dbfilename.substring(0, dbfilename.indexOf(".ws"));

        tmp = Character.toUpperCase(tmp.charAt(0))+tmp.substring(1);
        tmp = tmp.replace('_', ' ');

        if(!withdbsource)
            tmp = tmp.substring(0,tmp.indexOf(" [refseq"));

        return tmp;
    }

    /**
     * extract the database build date (modified date of .fasta file) from the database file name (format: Species_[refseq].mmddyyy.sdb)
     * e.g. input: rattus_norvegicus_[refseq].07_01_10.sdb
     * e.g. output: 07/01/10
     */
    public static String getDBBuild(String dbfilename){
        String tmp = dbfilename.substring(dbfilename.indexOf(".sdb")-8, dbfilename.indexOf(".sdb"));
        return tmp.replace("_", "/");
    }

    /**
     * extract the HID from locus_species
     * @param locus_species format is: HID|gene symbols_[species with spaces]
     * @return HID
     */
    public static String getHID(String locus_species){
        return locus_species.substring(0,locus_species.indexOf('|'));
    }

    /**
     * extract the locus from locus_species
     * @param locus_species format is: HID|gene symbols_[species with spaces]
     * @return HID|gene symbols
     */
    public static String getLocus(String locus_species){
        return locus_species.substring(0,locus_species.indexOf("_["));
    }

    /**
     * extract the species from locus_species
     * @param locus_species format is: HID|gene symbols_[species with spaces]
     * @return if nospace is true, return species with spaces replace by '_'
     *         otherwise, return species with spaces
     */
    public static String getSpecies(String locus_species, boolean nospace){
        String species = locus_species.substring(locus_species.indexOf('[')+1, locus_species.indexOf(']'));
        if(nospace)
            return species.replace(' ','_');
        else
            return species;
    }

    /**
     * build the following format: HID|gene symbol[species with or without spaces]
     * this format is used in by resultsFrame to display orthologous proteins with this to be used for column "HID|gene symbol[species]"
     * @return if nospace is true, return above format with species containing '_' in place of spaces
     *         otherwise, return above format with species containing spaces
     */
    public static String getHIDGeneSpecies(String locus_species, String genesymbol, boolean nospace){
        return getHID(locus_species)+"|"+genesymbol+"["+getSpecies(locus_species,nospace)+"]";
    }

    /**
     * build the following format: HID|species from locus_species
     * this format for the identifier of each protein sequence for ClustalW2
     * @param locus_species format is: HID|gene symbols_[species with spaces]
     * @return HID|species
     */
    public static String getHIDSpecies(String locus_species){
        return getHID(locus_species)+"|"+getSpecies(locus_species,true);
    }

    /**
     * from the database directory, find the homologene database (naming convention is homologene.buildnum.sdb)
     * @return the highest build number (buildnum)
     */
    public static int gethDBBuildNum() {
        File dbDir = new File(Analysis.dbFolderPath);
        int hbuild = 0;
        String[] dbs = dbDir.list();
        Arrays.sort(dbs);
        for (int i = 0; i < dbs.length; i++) {
            if (dbs[i].startsWith("homologene.") && dbs[i].endsWith(".sdb")) {
                int n = Integer.parseInt(dbs[i].substring(11, dbs[i].indexOf(".sdb")));
                if (n > hbuild) {
                    hbuild = n;
                }
            }
        }

        return hbuild;
    }

    /**
     * @return file name of the most current version of homologene (i.e. highest build number)
     */
    public static String gethDBfilename(){
        return "homologene."+Integer.toString(gethDBBuildNum())+".sdb";
    }

    /* process/thread handling methods */
    /**
     * executes a command in a separate process
     * @param commandStr command to be executed
     * @param separateThread if true, then start a new thread to execute the commands
     * @return true for successful execution; false otherwise
     */
    private static boolean sendToExec(final String[] commandStr) {
        try {
            //System.out.println(commandStr);
            Process p = Runtime.getRuntime().exec(commandStr);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String line;
            while((line = br.readLine()) != null){
                //System.out.println(line);
            }
            br.close();
            br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while((line = br.readLine()) != null){
                //System.out.println(line);
            }
            br.close();
            int exitVal = 0;
            try {
                exitVal = p.waitFor();
            } catch (InterruptedException e) {
            }

            if (exitVal == 0) {
                return true;
            }
        } catch (IOException e) {
            //e.printStackTrace();
        }

        return false;
    }

    /**
     * evaluates separateThread and pass on command value to main sendToExec method
     */
    public static boolean sendToExec(final String[] commandStr, final boolean separateThread) {

        if (separateThread) {
            new Thread() {
                @Override
                public void run() {
                    sendToExec(commandStr);
                }
            }.start();
        } else {
            return sendToExec(commandStr);
        }

        return true;
    }
    
    public static void submitAnalysisNoDisplay(String inputfilename, String outputdir, String soutputdir, int inputspeciesDB, int refspeciesDB){
        Utilities.updateSettings();
        System.out.println("Starting Cphosphosite with no display...");
        System.out.println("Input filename: "+inputfilename);
        System.out.println("Input species DB: "+inputspeciesDB);
        System.out.println("Input reference species DB: "+refspeciesDB);
        
        Analysis.alnOutFolderPath = outputdir;
        String input = Utilities.readFile(inputfilename);
        String speciespath = Analysis.db_species.get(inputspeciesDB);
        String dbpath = Analysis.db_target.get(refspeciesDB);
        Analysis analysis = new Analysis(speciespath, dbpath, null, false);

        String[] inputArray = input.split("\\s+");
        for (int i = 0; i < inputArray.length; i++)
            analysis.peptideInputs.add(inputArray[i].toUpperCase());

        analysis.analyzePeptides();

        //report running time as system out
        System.out.println(Utilities.getRuntimeStats(analysis.runningTime));
        
        //save outputs
        BufferedWriter inputFile = null;
        
        try {
        
            //phosphopeptide table
            inputFile = new BufferedWriter(new FileWriter(soutputdir+"phosphopeptide.txt"));
            inputFile.write("Phosphopeptide Sequence\tPhosphosite(s) ("+analysis.inputSpecies+")\n");

            for (int i = 0; i < analysis.peptideDB.size(); i++) {
                String s = analysis.peptideDB.get(i).get(0)+"\t";
                s += analysis.peptideDB.get(i).get(1)+"\n";
                inputFile.write(s);
            }
            inputFile.close();

            //phosphosite table
            inputFile = new BufferedWriter(new FileWriter(soutputdir+"phosphosite.txt"));
            inputFile.write("HID\tInput Species\tGene Symbol\tPhosphosite\n");
            for (int i = 0; i < analysis.phosphositeDB.size(); i++) {
                String s = Utilities.getHID(analysis.phosphositeDB.get(i).get(0))+"\t";
                s += Utilities.getSpecies(analysis.phosphositeDB.get(i).get(1),false)+"\t";
                s += analysis.phosphositeDB.get(i).get(3)+"\t";
                s += analysis.phosphositeDB.get(i).get(2)+"\n";
                inputFile.write(s);
            }
            inputFile.close();

            //phosphoprotein table
            inputFile = new BufferedWriter(new FileWriter(soutputdir+"phosphoprotein.txt"));
            inputFile.write("HID\tInput Species\tGene Symbol\tPhosphosite\n");
            for (int i = 0; i < analysis.protDB.size(); i++) {
                String s = Utilities.getHID(analysis.protDB.get(i).get(0))+"\t";
                s += Utilities.getSpecies(analysis.protDB.get(i).get(1),false)+"\t";
                s += analysis.protDB.get(i).get(2)+"\t";
                s += analysis.protDB.get(i).get(3)+"\t";
                s += analysis.protDB.get(i).get(4)+"\n";
                inputFile.write(s);
            }
            inputFile.close();

            //conservation score
            inputFile = new BufferedWriter(new FileWriter(soutputdir+"cscores.txt"));
            inputFile.write("HID\tGene Symbol ("+analysis.inputSpecies+")\tPhosphosite ("+analysis.inputSpecies+")\tSite Conservation Score\tMotif Conservation Score\tNumber of Orthologs\tCorresponding Phosphosites\n");
            HashMap<String, String> orthologsCount = new HashMap<String, String>(); //to be used in the next table
            for (int i = 0; i < analysis.locusList.size(); i++) {
                orthologsCount.put(analysis.locusList.get(i).get(0),analysis.locusList.get(i).get(1)); //for next table
            }
            for (int i = 0; i < analysis.cscoresDB.size(); i++) {
                String s = Utilities.getHID(analysis.cscoresDB.get(i).get(0))+"\t";
                s += analysis.phosphositeDB.get(i).get(3)+"\t";
                s += analysis.cscoresDB.get(i).get(1)+"\t";
                //round to three decimal places
                double a = Double.valueOf(analysis.cscoresDB.get(i).get(2));
                a = Math.round(a*1000.0)/1000.0;
                s += String.valueOf(a)+"\t";
                double b = Double.valueOf(analysis.cscoresDB.get(i).get(3));
                b = Math.round(b*1000.0)/1000.0;
                s += String.valueOf(b)+"\t";
                s += orthologsCount.get(analysis.cscoresDB.get(i).get(0))+"\t";
                s += analysis.phosphositeDB.get(i).get(4)+"\n";
                inputFile.write(s);
            }
            inputFile.close();
        } catch (IOException e){
            System.out.println("An error occured while trying to write to file.");
        }
    }
}
