package com.locaris.algorithm.beans;

public class TagPos {
	private String tagID;
	private float x;
	private float y;
	private float z;
	private float residual;
	private int panid;
	
	public float getResidual() {
		return residual;
	}
	public void setResidual(float residual) {
		this.residual = residual;
	}
	public float getX() {
		return x;
	}
	public void setX(float x) {
		this.x = x;
	}
	public float getY() {
		return y;
	}
	public void setY(float y) {
		this.y = y;
	}
	public float getZ() {
		return z;
	}
	public void setZ(float z) {
		this.z = z;
	}
	public int getPanid() {
		return panid;
	}
	public void setPanid(int panid) {
		this.panid = panid;
	}
	public String getTagID() {
		return tagID;
	}
	public void setTagID(String tagID) {
		this.tagID = tagID;
	}
	
	
	
}
