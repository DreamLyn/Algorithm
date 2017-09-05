package com.locaris.algorithm.location;

import com.locaris.algorithm.beans.Position;

/**
 * 接口定义了网格滤波算法，需实现此接口。
 * @author Mll
 */
public interface GridAlgorithmInterface {
	
	/**
	 * 使用网格初始化获取定位坐标及网格边界
	 * @param rdif------------>据差值
	 * @param devicePosition-->基站位置信息
	 * @param z_hat----------->高度差修正
	 * @return
	 */
	public Position pmfINITIALIZE(double[][] rdif,double[][] devicePosition,double z_hat);
	
	/**
	 * 使用网格滤波算法获取定位坐标及网格边界
	 * @param rdif-------------->距差值
	 * @param devicePosition---->基站位置信息
	 * @param z_hat------------->标签经验高度
	 * @param t----------------->最近两次定位的时间间隔t=（seq(i)-seq(i-1))*刷新频率
	 * @param position---------->上次定位输出结果
	 * @return
	 */
	public Position pmfTDOA(double[][] rdif,double[][] devicePosition,double z_hat,double t,Position positionOld);
}
