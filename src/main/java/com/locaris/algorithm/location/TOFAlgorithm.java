package com.locaris.algorithm.location;

import org.w3c.dom.ls.LSException;

import com.locaris.algorithm.beans.Position;

import Jama.Matrix;

public class TOFAlgorithm implements TOFAlgorithmInterface {
	private Position LSE(double[][] devicePosition, double[] range) {
		// TODO Auto-generated method stub
		int length = devicePosition.length;
		double[][] A = new double[devicePosition.length - 1][1];
		double[][] midb = new double[devicePosition.length - 1][1];
		double[][] B = new double[devicePosition.length - 1][1];
		for (int count = 0; count < length - 1; count++) {
			A[count][0] = 2 * (devicePosition[count + 1][1] - devicePosition[0][1]);
			midb[count][0] = range[0] * range[0] - range[count + 1] * range[count + 1];
			B[count][0] = midb[count][0] + devicePosition[count + 1][1] * devicePosition[count + 1][1]
					- devicePosition[0][1] * devicePosition[0][1];
		}
		Matrix aMatrix = new Matrix(A);
		Matrix bMatrix = new Matrix(B);
		Matrix yMatrix = aMatrix.transpose().times(aMatrix).inverse().times(aMatrix.transpose()).times(bMatrix);
		double y = yMatrix.get(0, 0);
		double x = Math.sqrt(range[0] * range[0] - (y - devicePosition[0][1]) * (y - devicePosition[0][1]));
		Position position = new Position();
		position.setX(x);
		position.setY(y);
		return position;
	}

	@Override
	public Position tofLSE(double[][] devicePosition, double[] range, double diff) {
		/**
		 * 对坐标进行旋转
		 */
		if (Math.abs(devicePosition[0][0] - devicePosition[1][0]) <= diff) {
			// 认为基站X轴相差很小,与Y轴平行.
			int length = devicePosition.length;
			double x=devicePosition[0][0];
			for (int count = 0; count < length; count++) {
				devicePosition[count][0] = 0;
			}
			Position position = this.LSE(devicePosition, range);
			position.setX(x);
			return position;
		} else {
			int length = devicePosition.length;
			// 斜率
			double slope = (devicePosition[0][1] - devicePosition[1][1])
					/ (devicePosition[0][0] - devicePosition[1][0]);
			double y0 = devicePosition[0][1] - devicePosition[0][0] * slope;
			// 基站到截距的距离
			double[] distance=new double[length];
			double[][] deviceP = new double[length][2];
			for (int count = 0; count < 2; count++) {
				distance[count]=Math.sqrt(devicePosition[count][0] * devicePosition[count][0]
						+ (devicePosition[count][1] - y0) * (devicePosition[count][1] - y0));
				deviceP[count][0]=0;
				deviceP[count][1] = distance[count]+y0;
			}
			Position position = this.LSE(deviceP, range);
			
			distance[0]=Math.abs(distance[0]-Math.abs(position.getY()-y0));
			distance[1]=Math.abs(distance[1]-Math.abs(position.getY()-y0));
			double x=(devicePosition[1][0]-devicePosition[0][0])*distance[0]/(distance[0]+distance[1])+devicePosition[0][0];
			double y=(devicePosition[1][1]-devicePosition[0][1])*distance[0]/(distance[0]+distance[1])+devicePosition[0][1];
			position.setX(x);
			position.setY(y);
			return position;
		}
	}

	
	
}
