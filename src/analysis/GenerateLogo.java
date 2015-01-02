package analysis;

import java.awt.*;
import javax.swing.*;
import java.util.*;

/**
 *
 * @author Boyang Zhao
 * @version 1.0
 * @description class for generating sequence logos
 */
public class GenerateLogo {

    private static class Logo extends JPanel {
        private final char character;
        private final double scale;
        private final Color charColor;
        public static Font charFont;
        public static int charheight;


        public Logo(char character, double scale, Color charColor) {
            this.character = character;
            this.scale = scale;
            this.charColor = charColor;
            setBorder(BorderFactory.createLineBorder(Color.blue));
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            //Graphics2D g2 = (Graphics2D) g;
            //java.awt.geom.AffineTransform at = new java.awt.geom.AffineTransform();
            //at.shear(0.2, 0.0);
            //g2.transform(at);
            //g2.setFont(new Font("Arial", Font.BOLD, getSize().width/25));
            //g2.setColor(Color.blue);
            //FontMetrics fm = this.getFontMetrics(charFont);
            //int height = fm.getAscent();
            //g2.drawString("A", 0, height);
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setFont(new Font("Arial", Font.BOLD, 100));
            g2.scale(0.5, scale);
            g2.setColor(charColor);
            g2.drawString(String.valueOf(character), 0, 72);
        }
    }

    /**
     * for display vertical label
     */
    private static class JVLabel extends JLabel {
        public JVLabel(String s) {
            super(s);
            setPreferredSize(new Dimension(30,100));
            setMinimumSize(new Dimension(30,100));
        }
        @Override
        public void paintComponent(Graphics g) {
            Graphics2D g2 = (Graphics2D) g;
            java.awt.geom.AffineTransform at = g2.getTransform();
            at.rotate(Math.toRadians(90));
            g2.setTransform(at);
            super.paintComponent(g);
            /*
            Graphics2D g2 = (Graphics2D) g;
            //g2.translate(30, 100);
            g2.rotate(Math.toRadians(-90));
            g2.drawString(this.getText(), 0, 0);*/
        }
    }

    public static void plot(char[][] seqAlign, java.awt.Window parent) {

        ConservationScore cs = new ConservationScore(seqAlign, true);

        //used for ordering from highest to lowest
        //find out what the highest sum is
        double highestSum = 0.0;
        for(int i = 0; i<cs.icsum.length; i++)
            highestSum = (cs.icsum[i] > highestSum) ? cs.icsum[i] : highestSum;

        //find out what the highest IC value is
        double highestIC = 0.0;
        for(int i = 0; i<cs.icmatrix.length; i++){
            for(int j = 0; j<cs.icmatrix[i].length; j++)
                highestIC = (cs.icmatrix[i][j] > highestIC) ? cs.icmatrix[i][j] : highestIC;
        }

        //metrics
        Logo.charheight = 72;
        Logo.charFont = new Font("Arial", Font.BOLD, 40);
        
        int width = 600;
        int height = 300;

        //create layouts
        JDialog seqlogo = new JDialog(parent);
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setSize(width, height);
        Font labelFont = new Font("Arial",Font.BOLD,12);

        //title
        JLabel titleLabel = new JLabel("Sequence Logo");
        titleLabel.setHorizontalAlignment(JLabel.CENTER);
        titleLabel.setFont(labelFont);
        mainPanel.add(titleLabel, BorderLayout.PAGE_START);

        //x-axis and y-axis
        JLabel xaxisLabel = new JLabel("residue positions");
        xaxisLabel.setHorizontalAlignment(JLabel.CENTER);
        xaxisLabel.setFont(labelFont);
        mainPanel.add(xaxisLabel, BorderLayout.PAGE_END);

        JVLabel yaxisLabel = new JVLabel("Bits");
        yaxisLabel.setVerticalAlignment(JLabel.CENTER);
        yaxisLabel.setFont(labelFont);
        mainPanel.add(yaxisLabel, BorderLayout.LINE_START);

        //plot
        JPanel plotPanel = new JPanel(new GridLayout(1,cs.icsum.length));
        for(int i = 0; i<cs.icsum.length; i++){
            JPanel poslogo = new JPanel();
            poslogo.setLayout(new BoxLayout(poslogo, BoxLayout.Y_AXIS));
            
            //put the IC values into a treemap so it's sorted (from largest to smallest)
            Comparator<Double> reversed = Collections.reverseOrder();
            TreeMap<Double, Character> ICval = new TreeMap<Double, Character>(reversed);
            for (int g = 0; g < cs.groupingsLabel.size(); g++)
                ICval.put(cs.icmatrix[i][g], cs.groupingsLabel.get(g));

            //add each of the character now
            Iterator iter = ICval.entrySet().iterator();
            while(iter.hasNext()){
                Map.Entry<Double, Character> entry = (Map.Entry<Double, Character>) iter.next();

                char c = entry.getValue();
                double scale = entry.getKey()/highestIC;
                
                poslogo.add(new Logo(c,1,cs.groupingsColor.get(c)));
            }
            

            poslogo.setBorder(BorderFactory.createLineBorder(Color.black));
            plotPanel.add(poslogo);
        }
        mainPanel.add(plotPanel, BorderLayout.CENTER);

        //add the main panel and display the dialog
        seqlogo.getContentPane().add(mainPanel);
        seqlogo.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        seqlogo.setLocationRelativeTo(parent);
        seqlogo.setSize(width,height);
        seqlogo.setVisible(true);

        /*
        //used for ordering from highest to lowest
        //find out what the highest sum is
        double highestSum = 0.0;
        for(int i = 0; i<cs.icsum.length; i++)
            highestSum = (cs.icsum[i] > highestSum) ? cs.icsum[i] : highestSum;

        //find out what the highest IC value is
        double highestIC = 0.0;
        for(int i = 0; i<cs.icmatrix.length; i++){
            for(int j = 0; j<cs.icmatrix[i].length; j++)
                highestIC = (cs.icmatrix[i][j] > highestIC) ? cs.icmatrix[i][j] : highestIC;
        }

        //metrics
        int maxheight_char = 130;
        int maxwidth_char = 40;
        Logo.charwidth = maxwidth_char;
        Logo.charheight = maxheight_char;
        Logo.charFont = new Font("Arial", Font.BOLD, 150);

        //scale - conservation from IC val to height
        int charscale = (int) (maxheight_char/highestIC);
        int maxheight_plot = (int) (charscale*highestSum);

        //generate logo
        
        for(int i = 0; i < cs.icmatrix.length; i++){
            int x = i*maxwidth_char;
            int y = 0;

            //put the IC values into a treemap so it's sorted (from largest to smallest)
            Comparator<Double> reversed = Collections.reverseOrder();
            TreeMap<Double, Character> ICval = new TreeMap<Double, Character>(reversed);
            for (int g = 0; g < cs.groupingsLabel.size(); g++)
                ICval.put(cs.icmatrix[i][g], cs.groupingsLabel.get(g));

            //now create the logo in same ordering as ICval
            //first: create filler spaces on the top for each column
            int height = maxheight_plot - (int) (charscale*cs.icsum[i]);
            JPanel panel = new JPanel(new GridLayout(1, 1));
            panel.add(new Logo(' ', 0, Color.WHITE));
            panel.setBounds(x, y, maxwidth_char, height);
            seqlogo.getContentPane().add(panel);
            y += height;
            
            //second: the characters
            Iterator iter = ICval.entrySet().iterator();
            while(iter.hasNext()){
                Map.Entry<Double, Character> entry = (Map.Entry<Double, Character>) iter.next();
                
                char c = entry.getValue();
                double scale = entry.getKey()/highestIC;
                height = (int) (scale*maxheight_char);
                
                panel = new JPanel(new GridLayout(1, 1));
                panel.add(new Logo(c,scale,cs.groupingsColor.get(c)));
                panel.setBounds(x, y, maxwidth_char, height);
                seqlogo.getContentPane().add(panel);
                y += height;
            }
        }*/
    }
}
