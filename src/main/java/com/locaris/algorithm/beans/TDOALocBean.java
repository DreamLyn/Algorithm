package com.locaris.algorithm.beans;

/**
 * TDOA数据接收
 * @author Lyn
 */
public class TDOALocBean implements Cloneable{
	private long deviceId;
	private long timestamp;

	public long getDeviceId() {
		return deviceId;
	}

	public void setDeviceId(long deviceId) {
		this.deviceId = deviceId;
	}

	public long getTimestamp() {
		return timestamp;
	}

	public void setTimestamp(long timestamp) {
		this.timestamp = timestamp;
	}
}
