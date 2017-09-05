package com.locaris.algorithm.beans;

public class DevicePos {
	private String deviceID;
	private float x;
	private float y;
	private float z;
	private byte updown;
	public String getDeviceID() {
		return deviceID;
	}
	public void setDeviceID(String deviceID) {
		this.deviceID = deviceID;
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
	public byte getUpdown() {
		return updown;
	}
	public void setUpdown(byte updown) {
		this.updown = updown;
	}

}
