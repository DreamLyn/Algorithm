package com.locaris.algorithm.beans;
/**
 * 类解释标签的xyz坐标
 * @author Lyn
 */
public class Position {
	private byte sequence;
	private double x;
	private double y;
	private double z;
	private double gridx_min;
	private double gridx_max;
	private double gridy_min;
	private double gridy_max;
	private int zoneid;
	public int getZoneid() {
		return zoneid;
	}
	public void setZoneid(int zoneid) {
		this.zoneid = zoneid;
	}
	public double getGridx_min() {
		return gridx_min;
	}
	public void setGridx_min(double gridx_min) {
		this.gridx_min = gridx_min;
	}
	public double getGridx_max() {
		return gridx_max;
	}
	public void setGridx_max(double gridx_max) {
		this.gridx_max = gridx_max;
	}
	public double getGridy_min() {
		return gridy_min;
	}
	public void setGridy_min(double gridy_min) {
		this.gridy_min = gridy_min;
	}
	public double getGridy_max() {
		return gridy_max;
	}
	public void setGridy_max(double gridy_max) {
		this.gridy_max = gridy_max;
	}
	public double getX() {
		return x;
	}
	public void setX(double x) {
		this.x = x;
	}
	public double getY() {
		return y;
	}
	public void setY(double y) {
		this.y = y;
	}
	public double getZ() {
		return z;
	}
	public void setZ(double z) {
		this.z = z;
	}
	public byte getSequence() {
		return sequence;
	}
	public void setSequence(byte sequence) {
		this.sequence = sequence;
	}
	
}
