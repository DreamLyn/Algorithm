package com.locaris.algorithm.residual;

import com.locaris.algorithm.beans.Position;

public class Residual {

    /**
     * 计算残差值
     *
     * @param rdif------------->距差
     * @param devicePosition--->基站位置信息
     * @param position--------->标签定位结果
     * @return
     */
    public static double calculateResidual(double[][] rdif, double[][] devicePosition, Position position) {
        double lengthToTag[] = new double[rdif.length];
        double length_diff[] = new double[rdif.length - 1];
        double allSum = 0;
        for (int count = 0; count < rdif.length; count++) {
            double doubleTemp = (devicePosition[count][0] - position.getX())
                    * (devicePosition[count][0] - position.getX())
                    + (devicePosition[count][1] - position.getY())
                    * (devicePosition[count][1] - position.getY());
            lengthToTag[count] = doubleTemp;
            if (count != 0) {
                length_diff[count - 1] = lengthToTag[count]
                        - lengthToTag[0];
                allSum += (length_diff[count - 1] - rdif[count - 1][0])
                        * (length_diff[count - 1] - rdif[count - 1][0]);
            }
        }
        return Math.sqrt(allSum / rdif.length);
    }
}
