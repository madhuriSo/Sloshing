package main;

import org.knowm.xchart.QuickChart;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.BitmapEncoder;

import java.awt.*;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

/**
 * Slosh visualization
 */
public class Visualizer {

    public static void main(String[] args) throws Exception {


        // Create chart
        final XYChart chart = QuickChart.getChart("Slosh Effect", "Tank ", "Height Of Liquid", "Water Level", new double[]{0}, new double[]{0});
        chart.getStyler().setChartBackgroundColor(Color.white);
        chart.getStyler().setPlotGridHorizontalLinesVisible(false);
        chart.getStyler().setPlotGridVerticalLinesVisible(false);

        //Set Max and min height

        chart.getStyler().setYAxisMax(0.04);
        chart.getStyler().setYAxisMin(-5.74E-03);

        final SwingWrapper<XYChart> sw = new SwingWrapper<XYChart>(chart);
        sw.displayChart();

     
         java.util.List<double[]> floatVal=readData();
        double[] currentval=new double[100];
        double[] location=new double[99];
        double q=1.0;

        for(int t=0;t<99;t++){
            location[t]=q;
            q++;
        }

        int i=2;
        while(i<20480){
            //Update chart
            currentval=floatVal.get(i);
            chart.updateXYSeries("Water Level",location,currentval,null);
            sw.repaintChart();
            if(i==15000){                       // Take screenshot at specified time
                BitmapEncoder.saveBitmapWithDPI(chart,"screenshot", BitmapEncoder.BitmapFormat.PNG,200);

            }

            i++;
            Thread.sleep(2);

        }

    }

    //Read data from the CSV file

    public static java.util.List<double[]> readData()
    {
        Path path = Paths.get("/Users/Madhuri/Documents/SLOSH/SMVisualizer/src/main/eta.csv");
        java.util.List<String[]> buf = new ArrayList<>();
        try {
            Files.lines(path).forEach(line -> {
                buf.add(line.split(","));
            });
        } catch(IOException e) {
            e.printStackTrace();
        }
        int size=buf.size();

        java.util.List<double[]> floatVal= new ArrayList<double[]>();


        for(int i=1;i<size;i++){
            double[] values= new double[99];
            for(int j=1;j<100;j++){
                values[j-1]=Double.parseDouble(buf.get(i)[j]);
            }
            floatVal.add(values);
        }
        return floatVal;

    }


    //Read data and multiply by bias to reduce size of graph

    public static java.util.List<double[]> readDatabias(Double bias)
    {
        Path path = Paths.get("/Users/Madhuri/Documents/SLOSH/SMVisualizer/src/main/eta.csv");
        java.util.List<String[]> buf = new ArrayList<>();
        try {
            Files.lines(path).forEach(line -> {
                buf.add(line.split(","));
            });
        } catch(IOException e) {
            e.printStackTrace();
        }
        int size=buf.size();

        java.util.List<double[]> floatVal= new ArrayList<double[]>();


        for(int i=1;i<size;i++){
            double[] values= new double[100];
            for(int j=1;j<100;j++){
                values[j-1]=(Double.parseDouble(buf.get(i)[j])*(bias));
            }

            floatVal.add(values);

        }
        return floatVal;

    }
}

