package gui;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.*;
import org.netbeans.swing.outline.*;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.table.*;
import java.awt.event.*;
import java.awt.Cursor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Arrays;

import analysis.*;

/*
 * resultsFrame.java
 *
 * Created on Jun 23, 2010, 9:54:17 AM
 */
/**
 *
 * @author Boyang Zhao
 */
public class resultsFrame extends javax.swing.JFrame {

    /* start of treetable model implementations for ortholog table */
    private class RenderData implements RenderDataProvider {

        private ArrayList<ArrayList<String> > locusList;
        
        public RenderData(final ArrayList<ArrayList<String> > locusIn){
            this.locusList = locusIn;
        }
        
        @Override
        public java.awt.Color getBackground(Object o){
            return null;
        }

        @Override
        public String getDisplayName(Object o){
            Integer[] entryIndex = (Integer[]) o;

            if(entryIndex[0] == 1){
                return locusList.get(entryIndex[1]).get(0)+" ("+locusList.get(entryIndex[1]).get(1)+")";
            } else
                return "";
        }

        @Override
        public java.awt.Color getForeground(Object o){
            return null;
        }

        @Override
        public javax.swing.Icon getIcon(Object o){
            return null;
        }

        @Override
        public String getTooltipText(Object o){
            return null;
        }

        @Override
        public boolean isHtmlDisplayName(Object o){
            return false;
        }
    }
    
    private static class OrthologTreeModel implements TreeModel {

        //entryIndex data structure: Array of Integers
        //[0] represents the level; [1] represents the index of entry in database
        //this structure is made to optimize TreeTable display, so don't have to parse
        //through the orthologDB

        //for entryIndex[0]=0: at root level;
        //for entryIndex[0]=1: at locus level; child is for the list of locus
        //                     entryIndex[1] represents index for dbMain
        //for entryIndex[0]=2: at orthologs level; child is the list of orthologs with same locus
        //                     entryIndex[1] represents index for locusList

        private final ArrayList<ArrayList<String> > dbMain;
        private final ArrayList<ArrayList<String> > locusList;

        public OrthologTreeModel(final ArrayList<ArrayList<String> > dbIn, final ArrayList<ArrayList<String> > locusIn){
            this.dbMain = dbIn;
            this.locusList = locusIn;
        }

        @Override
        public void addTreeModelListener(javax.swing.event.TreeModelListener lstnr){
        }

        @Override
        public Object getChild(Object parent, int index){
            Integer[] par = (Integer[]) parent;
            Integer[] entryIndex = new Integer[2];
            
            if(par[0] == 0){//at the root level
                entryIndex[0] = 1; //child of root is list of locus
                entryIndex[1] = index; //index of locusList
            } else if(par[0] == 1){//at the locus listing level
                entryIndex[0] = 2; //child of locus is list of orthologs
                entryIndex[1] = Integer.parseInt(locusList.get(par[1]).get(2))+index; //index of dbMain
            }

            return entryIndex;
        }

        @Override
        public int getChildCount(Object parent){

            Integer[] entryIndex = (Integer[]) parent;
            if(entryIndex[0] == 0)
                return locusList.size();
            else if(entryIndex[0] == 1)
                return Integer.parseInt(locusList.get(entryIndex[1]).get(1));
            else
                return 0;
        }

        @Override
        public int getIndexOfChild(Object parent, Object child){
            Integer[] entryIndexP = (Integer[]) parent;
            Integer[] entryIndexC = (Integer[]) child;
            
            if(entryIndexC[0] == 1)
                return entryIndexC[1];
            else if(entryIndexC[0] == 2)
                return entryIndexC[1]-Integer.parseInt(locusList.get(entryIndexP[1]).get(2));

            return -1;
        }

        @Override
        public Object getRoot(){
            Integer[] entryIndex = {0,0};
            return entryIndex;
        }

        @Override
        public boolean isLeaf(Object node){
            Integer[] entryIndex = (Integer[]) node;
            if(entryIndex[0] == 2)
                return true;
            return false;
        }

        @Override
        public void removeTreeModelListener(javax.swing.event.TreeModelListener lstnr){
            //do nothing
        }

        @Override
        public void valueForPathChanged(javax.swing.tree.TreePath path, Object newValue){
            //do nothing
        }
    }

    private class OrthologRowModel implements RowModel {

        private final ArrayList<ArrayList<String> > dbMain;
        private final ArrayList<ArrayList<String> > locusList;

        public OrthologRowModel(final ArrayList<ArrayList<String>> dbIn, final ArrayList<ArrayList<String>> locusIn){
            this.dbMain = dbIn;
            this.locusList = locusIn;
        }

        @Override
        public Class getColumnClass(int column){
            return String.class;
        }

        @Override
        public int getColumnCount() {
            return 4;
        }

        @Override
        public String getColumnName(int column){
            switch (column) {
                case 0:
                    return "HID|gene symbol[species]";
                case 1:
                    return "RefSeq Accession Number";
                case 2:
                    return "Gene Symbol";
                case 3:
                    return "Protein Sequence";
            }

            return null;
        }

        @Override
        public Object getValueFor(Object node, int column) {
            Integer[] entryIndex = (Integer[]) node;

            if (entryIndex[0] == 1) {
                return "";
            } else if (entryIndex[0] == 2) {
                if (column == 0) {
                    return Utilities.getHIDGeneSpecies(dbMain.get(entryIndex[1]).get(1), dbMain.get(entryIndex[1]).get(3), false);
                } else {
                    return dbMain.get(entryIndex[1]).get(column + 1);
                }
            }

            return null;
        }

        @Override
        public boolean isCellEditable(Object node, int column) {
            return false;
        }

        @Override
        public void setValueFor(Object node, int column, Object value) {
            //do nothing for now
        }
    }
    /* end of treetable model implementations for ortholog table */

    /* start of jbutton in jtable classes */
    private class JTblBtnRenderer extends JButton implements TableCellRenderer {

        public JTblBtnRenderer(String text) {
            this.setText(text);
        }

        @Override
        public java.awt.Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
            return this;
        }
    }

    private class JTblBtnMouseListener implements MouseListener {

        private JTable table;

        public JTblBtnMouseListener(JTable table_i){
            table = table_i;
        }
        
        @Override
        public void mouseClicked(MouseEvent e) {
        }

        @Override
        public void mouseExited(MouseEvent e) {
            setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        }

        @Override
        public void mouseEntered(MouseEvent e) {
            int rowIndex = table.rowAtPoint(e.getPoint());
            int colIndex = table.columnAtPoint(e.getPoint());

            if (!(table.getValueAt(rowIndex, colIndex) instanceof JButton)) {
                setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
            } else {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            }
        }

        @Override
        public void mouseReleased(MouseEvent e) {
        }

        @Override
        public void mousePressed(MouseEvent e) {
        }
    }
    /* end of jbutton in jtable classes */

    /* constructor */
    public resultsFrame(final Analysis analysisI) {
        this.analysis = analysisI;
        this.resultsFrame = this;
        populateTables();
        initComponents();
        postinit();
    }

    private void postinit() {
        //for seq align table
        //execute external programs to display seq alignment or tree when cell is clicked
        seqalignTable.getColumn("View Sequence Alignment").setCellRenderer(new JTblBtnRenderer("view seq alignment"));
        seqalignTable.getColumn("View Phylogenetic Tree").setCellRenderer(new JTblBtnRenderer("view tree"));
        seqalignTable.addMouseListener(new JTblBtnMouseListener(seqalignTable) {

            @Override
            public void mouseClicked(MouseEvent e) {
                int rowIndex = seqalignTable.rowAtPoint(e.getPoint());
                int colIndex = seqalignTable.columnAtPoint(e.getPoint());
                String locus = analysis.locusList.get(rowIndex).get(0);
                String HID = Utilities.getHID(locus);
                
                switch(colIndex){
                    case 2:
                        //view seq alignment
                        String[] ClustalCommand = {Analysis.exeFolderPath+"clustalx", "-infile="+Analysis.alnOutFolderPath+HID+".aln"};
                        
                        Utilities.sendToExec(ClustalCommand, true);
                        break;
                    case 3:
                        //view tree
                        String[] NJPlotCommand = {Analysis.exeFolderPath+"njplot", Analysis.alnOutFolderPath+HID+".ph"};
                        
                        Utilities.sendToExec(NJPlotCommand, true);
                        break;
                }
            }
        });

        //for conservation score table
        /* view logo feature; not yet implemented
        cscoresTable.getColumn("View Sequence Logo").setCellRenderer(new JTblBtnRenderer("view logo"));
        cscoresTable.addMouseListener(new JTblBtnMouseListener(cscoresTable) {

            @Override
            public void mouseClicked(MouseEvent e) {
                int rowIndex = cscoresTable.rowAtPoint(e.getPoint());
                int colIndex = cscoresTable.columnAtPoint(e.getPoint());

                if (colIndex == 7) {
                    GenerateLogo.plot(analysis.phosphositeAlignDB.get(rowIndex), resultsFrame);
                }
            }
        });*/

        //sets the width of tables
        int w;
        w = overviewTable.getWidth();
        overviewTable.getColumnModel().getColumn(0).setPreferredWidth(70*w/100);
        overviewTable.getColumnModel().getColumn(1).setPreferredWidth(30*w/100);

        w = phosphoproteinTable.getWidth();
        phosphoproteinTable.getColumnModel().getColumn(0).setPreferredWidth(10*w/100);
        phosphoproteinTable.getColumnModel().getColumn(1).setPreferredWidth(20*w/100);
        phosphoproteinTable.getColumnModel().getColumn(2).setPreferredWidth(25*w/100);
        phosphoproteinTable.getColumnModel().getColumn(3).setPreferredWidth(15*w/100);
        phosphoproteinTable.getColumnModel().getColumn(4).setPreferredWidth(30*w/100);

        w = orthologTreeTable.getWidth();
        orthologTreeTable.getColumnModel().getColumn(0).setPreferredWidth(35*w/100);
        orthologTreeTable.getColumnModel().getColumn(1).setPreferredWidth(30*w/100);
        orthologTreeTable.getColumnModel().getColumn(2).setPreferredWidth(15*w/100);
        orthologTreeTable.getColumnModel().getColumn(3).setPreferredWidth(10*w/100);
        orthologTreeTable.getColumnModel().getColumn(4).setPreferredWidth(10*w/100);
    }

    private void populateTables() {
        overviewTableRow = new Object[][]{
                    {"input species database", Utilities.getDBName(Analysis.protFileName, true)},
                    {"database to be searched against", Utilities.getDBName(Analysis.orthoFileName, true)},
                    {"number of input phosphopeptides", analysis.peptideDB.size()},
                    {"number of phosphopeptides with ambigious protein and/or site identification", analysis.peptideDBambi.size()},
                    {"number of unique phosphosites", analysis.phosphositeDB.size()},
                    {"number of unique phosphoproteins", analysis.protDB.size()}
                };

        //phosphopeptide table
        phosphopeptideRow = new Object[analysis.peptideDB.size()][2];
        for (int i = 0; i < analysis.peptideDB.size(); i++) {
            phosphopeptideRow[i][0] = analysis.peptideDB.get(i).get(0);
            phosphopeptideRow[i][1] = analysis.peptideDB.get(i).get(1);
        }
        phosphopeptideCol = new String[]{"Phosphopeptide Sequence","Phosphosite(s) ("+analysis.inputSpecies+")"};

        //phosphosite table
        phosphositeRow = new Object[analysis.phosphositeDB.size()][4];
        for (int i = 0; i < analysis.phosphositeDB.size(); i++) {
            phosphositeRow[i][0] = Utilities.getHID(analysis.phosphositeDB.get(i).get(0));
            phosphositeRow[i][1] = Utilities.getSpecies(analysis.phosphositeDB.get(i).get(1),false);
            phosphositeRow[i][2] = analysis.phosphositeDB.get(i).get(3);
            phosphositeRow[i][3] = analysis.phosphositeDB.get(i).get(2);
        }
        phosphositeCol = new String[]{"HID","Input Species","Gene Symbol","Phosphosite"};
        
        //phosphoprotein table
        proteinDBRow = new Object[analysis.protDB.size()][5];
        for (int i = 0; i < analysis.protDB.size(); i++) {
            proteinDBRow[i][0] = Utilities.getHID(analysis.protDB.get(i).get(0));
            proteinDBRow[i][1] = Utilities.getSpecies(analysis.protDB.get(i).get(1),false);
            proteinDBRow[i][2] = analysis.protDB.get(i).get(2);
            proteinDBRow[i][3] = analysis.protDB.get(i).get(3);
            proteinDBRow[i][4] = analysis.protDB.get(i).get(4);
        }
        proteinDBCol = new String[]{"HID","Input Species","RefSeq Accession Number", "Gene Symbol", "Protein Sequence"};

        //ortholog table
        orthologMdl = DefaultOutlineModel.createOutlineModel(
                new OrthologTreeModel(analysis.orthologDB, analysis.locusList), new OrthologRowModel(analysis.orthologDB, analysis.locusList), true, "HID|all gene symbols");
        orthologRenderData = new RenderData(analysis.locusList);
        
        //sequence alignment
        seqalignRow = new Object[analysis.locusList.size()][4];
        HashMap<String, String> orthologsCount = new HashMap<String, String>(); //to be used in the next table
        for (int i = 0; i < analysis.locusList.size(); i++) {
            seqalignRow[i][0] = Utilities.getHID(analysis.locusList.get(i).get(0));
            seqalignRow[i][1] = analysis.orthologDB.get(Integer.parseInt(analysis.locusList.get(i).get(2))).get(3);
            seqalignRow[i][2] = new JButton();
            seqalignRow[i][3] = new JButton();
            orthologsCount.put(analysis.locusList.get(i).get(0),analysis.locusList.get(i).get(1)); //for next table
        }
        seqalignCol = new String[]{"HID","Gene Symbol ("+analysis.inputSpecies+")","View Sequence Alignment","View Phylogenetic Tree"};

        //conservation score
        cscoreRow = new Object[analysis.cscoresDB.size()][9];
        for (int i = 0; i < analysis.cscoresDB.size(); i++) {
            cscoreRow[i][0] = Utilities.getHID(analysis.cscoresDB.get(i).get(0));
            cscoreRow[i][1] = analysis.phosphositeDB.get(i).get(3);
            cscoreRow[i][2] = analysis.cscoresDB.get(i).get(1);
            //round to three decimal places
            double a = Double.valueOf(analysis.cscoresDB.get(i).get(2));
            a = Math.round(a*1e3)/1e3;
            cscoreRow[i][3] = String.valueOf(a);
            double b = Double.valueOf(analysis.cscoresDB.get(i).get(3));
            b = Math.round(b*1e3)/1e3;
            cscoreRow[i][4] = String.valueOf(b);
            cscoreRow[i][5] = orthologsCount.get(analysis.cscoresDB.get(i).get(0));
            cscoreRow[i][6] = analysis.phosphositeDB.get(i).get(4);
            cscoreRow[i][7] = analysis.cscoresDB.get(i).get(4);
            String hashkey = analysis.phosphositeDB.get(i).subList(0,4).toString().replaceAll("(^\\[|\\]$)","").replace(", ","\t");
            cscoreRow[i][8] = analysis.phosphositeSList.get(hashkey);
        }
        cscoreCol = new String[]{"HID","Gene Symbol ("+analysis.inputSpecies+")","Phosphosite ("+analysis.inputSpecies+")",
            "Site Conservation Score","Motif Conservation Score","Number of Orthologs","Corresponding Phosphosites","Conservation Scores Breakdown",
            "Input Indices"};
    }
    
    private void saveTable(Object[][] table, String[] tableCol, String filename, Boolean... options) {
        //export table and save at provided path
        String inputFilePath = Utilities.saveFile(this, filename);
        FileWriter f = null;
        Boolean outFormatCollapsed = true;
        if(options.length > 0){
            //first option: if boolean is false, then will output non-unique rows for each value in the last column of the table
            //for example, if the last column value is 2,4,6; then instead of writing the row once, will write the row three times with different
            //last column values
            //this was used for cscores_mapping.csv
            if (options[0] instanceof Boolean) {
                outFormatCollapsed = options[0];
            }
        }
        
        try {
            f = new FileWriter(inputFilePath);
            
            //write table header
            if(tableCol.length > 0){
                String row = "";
                for(String colname : tableCol){
                    row += colname+",";
                }
                row = row.substring(0,row.length());
                f.write(row+"\n");
            }
            
            //write table values
            for(int i = 0; i < table.length; i++){
                
                String[] lastColValList = {table[i][table[i].length-1].toString()};
                if(!outFormatCollapsed){
                    //will output expanded version
                    lastColValList = lastColValList[0].split(";");
                }
                
                String rowToRepeat = "";
                String row = "";
                for(int j = 0; j < table[i].length; j++){
                    if(!(table[i][j] instanceof JComponent)){
                        if(j != 0)
                            row += ",";

                        if(j == table[i].length-1){
                            if(lastColValList.length > 1){
                                rowToRepeat = row;
                                row += lastColValList[0];
                            } else {
                                row += table[i][j].toString();
                            }
                        } else {
                            row += table[i][j].toString();
                        }
                    }
                }
                f.write(row+"\n");
                
                if(lastColValList.length > 1){
                    for(int j = 1; j < lastColValList.length; j++){
                        row = rowToRepeat+lastColValList[j];
                        f.write(row+"\n");
                    }
                }
            }
        } catch (IOException e) {
        } finally {
            try {
                if (f != null)
                    f.close();
            } catch(IOException e){}
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        resultsJPanel = new javax.swing.JTabbedPane();
        overviewPanel = new javax.swing.JPanel();
        jScrollPane4 = new javax.swing.JScrollPane();
        overviewTable = new javax.swing.JTable();
        phosphopeptidePanel = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        phosphopeptideTable = new javax.swing.JTable();
        phosphositePanel = new javax.swing.JPanel();
        jScrollPane5 = new javax.swing.JScrollPane();
        phosphositeTable = new javax.swing.JTable();
        phosphoproteinPanel = new javax.swing.JPanel();
        jScrollPane2 = new javax.swing.JScrollPane();
        phosphoproteinTable = new javax.swing.JTable();
        orthologPanel = new javax.swing.JPanel();
        jScrollPane3 = new javax.swing.JScrollPane();
        orthologTreeTable = new org.netbeans.swing.outline.Outline();
        seqalignPanel = new javax.swing.JPanel();
        jScrollPane6 = new javax.swing.JScrollPane();
        seqalignTable = new javax.swing.JTable();
        cscoresPanel = new javax.swing.JPanel();
        jScrollPane7 = new javax.swing.JScrollPane();
        cscoresTable = new javax.swing.JTable();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        exportMenuItem = new javax.swing.JMenuItem();
        jSeparator1 = new javax.swing.JPopupMenu.Separator();
        closeMenuItem = new javax.swing.JMenuItem();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Results - CPhos");
        setIconImage(new ImageIcon(getClass().getResource("/resources/NHLBI.jpg")).getImage());

        resultsJPanel.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                resultsJPanelStateChanged(evt);
            }
        });

        overviewTable.setModel(new javax.swing.table.DefaultTableModel(
            overviewTableRow,
            new String [] {
                "Parameter", "Value"
            }
        ));
        jScrollPane4.setViewportView(overviewTable);

        javax.swing.GroupLayout overviewPanelLayout = new javax.swing.GroupLayout(overviewPanel);
        overviewPanel.setLayout(overviewPanelLayout);
        overviewPanelLayout.setHorizontalGroup(
            overviewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane4, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
        );
        overviewPanelLayout.setVerticalGroup(
            overviewPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane4, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
        );

        resultsJPanel.addTab("Overview", overviewPanel);

        phosphopeptideTable.setModel(new javax.swing.table.DefaultTableModel(
            phosphopeptideRow,
            phosphopeptideCol
        ){
            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return false;
            }
        }
    );
    jScrollPane1.setViewportView(phosphopeptideTable);

    javax.swing.GroupLayout phosphopeptidePanelLayout = new javax.swing.GroupLayout(phosphopeptidePanel);
    phosphopeptidePanel.setLayout(phosphopeptidePanelLayout);
    phosphopeptidePanelLayout.setHorizontalGroup(
        phosphopeptidePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    phosphopeptidePanelLayout.setVerticalGroup(
        phosphopeptidePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Phosphopeptides", phosphopeptidePanel);

    phosphositeTable.setModel(new javax.swing.table.DefaultTableModel(
        phosphositeRow,
        phosphositeCol
    ){
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }
    }
    );
    jScrollPane5.setViewportView(phosphositeTable);

    javax.swing.GroupLayout phosphositePanelLayout = new javax.swing.GroupLayout(phosphositePanel);
    phosphositePanel.setLayout(phosphositePanelLayout);
    phosphositePanelLayout.setHorizontalGroup(
        phosphositePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane5, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    phosphositePanelLayout.setVerticalGroup(
        phosphositePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane5, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Phosphosites", phosphositePanel);

    phosphoproteinTable.setModel(new javax.swing.table.DefaultTableModel(
        proteinDBRow,
        proteinDBCol
    ){
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }
    }
    );
    jScrollPane2.setViewportView(phosphoproteinTable);

    javax.swing.GroupLayout phosphoproteinPanelLayout = new javax.swing.GroupLayout(phosphoproteinPanel);
    phosphoproteinPanel.setLayout(phosphoproteinPanelLayout);
    phosphoproteinPanelLayout.setHorizontalGroup(
        phosphoproteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    phosphoproteinPanelLayout.setVerticalGroup(
        phosphoproteinPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane2, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Phosphoproteins", phosphoproteinPanel);

    orthologTreeTable.setRenderDataProvider(orthologRenderData);
    orthologTreeTable.setRootVisible(false);
    orthologTreeTable.setModel(orthologMdl);
    jScrollPane3.setViewportView(orthologTreeTable);

    javax.swing.GroupLayout orthologPanelLayout = new javax.swing.GroupLayout(orthologPanel);
    orthologPanel.setLayout(orthologPanelLayout);
    orthologPanelLayout.setHorizontalGroup(
        orthologPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    orthologPanelLayout.setVerticalGroup(
        orthologPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane3, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Orthologs", orthologPanel);

    seqalignTable.setModel(new javax.swing.table.DefaultTableModel(
        seqalignRow,
        seqalignCol
    ){
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }
    }
    );
    jScrollPane6.setViewportView(seqalignTable);

    javax.swing.GroupLayout seqalignPanelLayout = new javax.swing.GroupLayout(seqalignPanel);
    seqalignPanel.setLayout(seqalignPanelLayout);
    seqalignPanelLayout.setHorizontalGroup(
        seqalignPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    seqalignPanelLayout.setVerticalGroup(
        seqalignPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Sequence Alignments", seqalignPanel);

    cscoresTable.setModel(new javax.swing.table.DefaultTableModel(
        cscoreRow,
        cscoreCol
    ){
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return false;
        }
    }
    );
    jScrollPane7.setViewportView(cscoresTable);

    javax.swing.GroupLayout cscoresPanelLayout = new javax.swing.GroupLayout(cscoresPanel);
    cscoresPanel.setLayout(cscoresPanelLayout);
    cscoresPanelLayout.setHorizontalGroup(
        cscoresPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane7, javax.swing.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );
    cscoresPanelLayout.setVerticalGroup(
        cscoresPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addComponent(jScrollPane7, javax.swing.GroupLayout.DEFAULT_SIZE, 423, Short.MAX_VALUE)
    );

    resultsJPanel.addTab("Conservation Scores", cscoresPanel);

    jMenu1.setText("File");

    exportMenuItem.setText("Export current tab as .csv");
    exportMenuItem.addActionListener(new java.awt.event.ActionListener() {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            exportMenuItemActionPerformed(evt);
        }
    });
    jMenu1.add(exportMenuItem);
    jMenu1.add(jSeparator1);

    closeMenuItem.setText("Close");
    closeMenuItem.addActionListener(new java.awt.event.ActionListener() {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            closeMenuItemActionPerformed(evt);
        }
    });
    jMenu1.add(closeMenuItem);

    jMenuBar1.add(jMenu1);

    setJMenuBar(jMenuBar1);

    javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
    getContentPane().setLayout(layout);
    layout.setHorizontalGroup(
        layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addGroup(layout.createSequentialGroup()
            .addContainerGap()
            .addComponent(resultsJPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 858, Short.MAX_VALUE)
            .addContainerGap())
    );
    layout.setVerticalGroup(
        layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
        .addGroup(layout.createSequentialGroup()
            .addContainerGap()
            .addComponent(resultsJPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 451, Short.MAX_VALUE)
            .addContainerGap())
    );

    pack();
    }// </editor-fold>//GEN-END:initComponents

    private void closeMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_closeMenuItemActionPerformed
        dispose();
    }//GEN-LAST:event_closeMenuItemActionPerformed

    private void exportMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportMenuItemActionPerformed
        //export data as Excel spreadsheets
        int ctab = resultsJPanel.getSelectedIndex();
        switch(ctab){
            case 0:
                saveTable(overviewTableRow, new String[0], "./overviews.csv");
                break;
            case 1:
                saveTable(phosphopeptideRow, phosphopeptideCol, "./phosphopeptides.csv");
                break;
            case 2:
                saveTable(phosphositeRow, phosphositeCol, "./phosphosites.csv");
                break;
            case 3:
                saveTable(proteinDBRow, proteinDBCol, "./phosphoproteins.csv");
                break;
            case 4:
                //orthologs
                break;
            case 5:
                //sequence alignment
                break;
            case 6:
                saveTable(cscoreRow, cscoreCol, "./cscores.csv");
                saveTable(cscoreRow, cscoreCol, "./cscores_mapping.csv", false);
                break;
        }
    }//GEN-LAST:event_exportMenuItemActionPerformed

    private void resultsJPanelStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_resultsJPanelStateChanged
        int ctab = resultsJPanel.getSelectedIndex();
        switch(ctab){
            case 0:
            case 1:
            case 2:
            case 3:
            case 6:
                exportMenuItem.setEnabled(true);
                break;
            case 4:
            case 5:
                exportMenuItem.setEnabled(false);
                break;
        }
    }//GEN-LAST:event_resultsJPanelStateChanged
    /**
     * @param args the command line arguments
     */
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem closeMenuItem;
    private javax.swing.JPanel cscoresPanel;
    private javax.swing.JTable cscoresTable;
    private javax.swing.JMenuItem exportMenuItem;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JScrollPane jScrollPane4;
    private javax.swing.JScrollPane jScrollPane5;
    private javax.swing.JScrollPane jScrollPane6;
    private javax.swing.JScrollPane jScrollPane7;
    private javax.swing.JPopupMenu.Separator jSeparator1;
    private javax.swing.JPanel orthologPanel;
    private org.netbeans.swing.outline.Outline orthologTreeTable;
    private javax.swing.JPanel overviewPanel;
    private javax.swing.JTable overviewTable;
    private javax.swing.JPanel phosphopeptidePanel;
    private javax.swing.JTable phosphopeptideTable;
    private javax.swing.JPanel phosphoproteinPanel;
    private javax.swing.JTable phosphoproteinTable;
    private javax.swing.JPanel phosphositePanel;
    private javax.swing.JTable phosphositeTable;
    private javax.swing.JTabbedPane resultsJPanel;
    private javax.swing.JPanel seqalignPanel;
    private javax.swing.JTable seqalignTable;
    // End of variables declaration//GEN-END:variables
    private Object[][] phosphopeptideRow;
    private String[] phosphopeptideCol;
    private Object[][] phosphositeRow;
    private String[] phosphositeCol;
    private Object[][] proteinDBRow;
    private String[] proteinDBCol;
    private Object[][] seqalignRow;
    private String[] seqalignCol;
    private Object[][] cscoreRow;
    private String[] cscoreCol;
    private OutlineModel orthologMdl;
    private RenderData orthologRenderData;
    private Object[][] overviewTableRow;
    private final Analysis analysis;
    private JFrame resultsFrame;
}
