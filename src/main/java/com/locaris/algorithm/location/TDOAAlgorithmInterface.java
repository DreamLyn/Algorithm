package com.locaris.algorithm.location;

import com.locaris.algorithm.beans.Position;

/**
 * TDOA定位算法接口，在此接口所含算法中，所有距离信息单位都为米（m）
 */
public interface TDOAAlgorithmInterface {
	final long Max_TS = 1099511627776L;
	/**
	 * 基于chen的TDOA定位算法。
	 * 作者：李亚楠 张绘军
	 * @param rdif:据差值
	 * @param devicePosition:基站位置信息
	 *                      四列分别为 X坐标，Y坐标，Z坐标，平方和
	 * @param z_hat:高度差修正
	 * @return
	 */
	public Position chenTDOA(double[][] rdif, double[][] devicePosition, double z_hat);

	/**
	 * 基于srd的TDOA定位算法
	 * 作者：张绘军
	 * @param rdif
	 * @param devicePosition
	 * @param z_hat
	 * @return
	 */
	public Position srdTDOA(double[][] rdif, double[][] devicePosition, double z_hat);

	/**
	 * 用于狭长（long and narrow）环境的两基站一维TDOA定位
	 * 作者：张绘军
	 * @param rdif
	 * @param devicePosition
	 * @param z_hat
	 * @return
	 */
	public Position LnOneDim(double rdif, double[][] devicePosition, double z_hat);
}