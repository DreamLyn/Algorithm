package com.locaris.algorithm.location;

import com.locaris.algorithm.beans.Position;

public interface TOFAlgorithmInterface {
	/**
	 * TOF一维定位算法
	 * @param devicePosition
	 * @param range
	 * @param diff-->X轴相差多少认为X轴一样.
	 * @return
	 */
	public Position tofLSE(double[][] devicePosition,double[] range,double diff);
	
}
