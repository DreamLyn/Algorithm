package com.locaris.algorithm.location;

import com.locaris.algorithm.beans.Position;

/**
 * 接口定义了TDOA定位算法，需实现此接口。
 * 
 * @author Lyn
 */
public interface TDOAAlgorithmInterface {
	final long Max_TS = 1099511627776L;

	/**
	 * 使用chen算法获取定位坐标
	 * 
	 * @param rdif------------>据差值
	 * @param devicePosition-->基站位置信息
	 * @param z_hat----------->高度差修正
	 * @return
	 */
	public Position chenTDOA(double[][] rdif, double[][] devicePosition, double z_hat);

	/**
	 * 使用一维定位算法
	 * 
	 * @param rdif
	 * @param devicePosition
	 * @param z_hat
	 * @return
	 */
	public Position srdTDOA(double[][] rdif, double[][] devicePosition, double z_hat);

	public Position oneDim(double rdif, double[][] devicePosition, double z_hat);

	public Position LnOneDim(double rdif, double[][] devicePosition, double z_hat);
}
