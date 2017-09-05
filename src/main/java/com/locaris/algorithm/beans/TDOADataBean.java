package com.locaris.algorithm.beans;

/**
 * TDOA数据接收
 * @author Lyn
 */
public class TDOADataBean implements Cloneable{
	private byte dataIndex;
	private long deviceID;
	private long tagId;
	private byte locIndex;
	private long timeStamp;
	private byte syncIndex;
	private byte tagType;
	private short tagState;
	private short slotIndex;
	private boolean lastUse=false;
	private byte updown;
	public byte getUpdown() {
		return updown;
	}
	public void setUpdown(byte updown) {
		this.updown = updown;
	}
	public boolean isLastUse() {
		return lastUse;
	}
	public void setLastUse(boolean lastUse) {
		this.lastUse = lastUse;
	}
	public short getSlotIndex() {
		return slotIndex;
	}
	public void setSlotIndex(short slotIndex) {
		this.slotIndex = slotIndex;
	}
	public byte getTagType() {
		return tagType;
	}
	public void setTagType(byte tagType) {
		this.tagType = tagType;
	}
	private byte status;
	private float firstPathRssi;
	private float rssi;
	/**
	 * 上面为接收到的定位数据,下面为引擎所需数据.
	 */
	private float x;
	private float y;
	private float z;

	public byte getDataIndex() {
		return dataIndex;
	}
	public void setDataIndex(byte dataIndex) {
		this.dataIndex = dataIndex;
	}
	
	public short getTagState() {
		return tagState;
	}
	public void setTagState(short tagState) {
		this.tagState = tagState;
	}
	public byte getLocIndex() {
		return locIndex;
	}
	public void setLocIndex(byte locIndex) {
		this.locIndex = locIndex;
	}
	public byte getSyncIndex() {
		return syncIndex;
	}
	public void setSyncIndex(byte syncIndex) {
		this.syncIndex = syncIndex;
	}
	@Override  
    public Object clone() throws CloneNotSupportedException  
    {  
        return super.clone();  
    }  
	public long getDeviceID() {
		return deviceID;
	}
	public void setDeviceID(long deviceID) {
		this.deviceID = deviceID;
	}
	public long getTagId() {
		return tagId;
	}
	public void setTagId(long tagId) {
		this.tagId = tagId;
	}
	public long getTimeStamp() {
		return timeStamp;
	}
	public void setTimeStamp(long timeStamp) {
		this.timeStamp = timeStamp;
	}
	public byte getStatus() {
		return status;
	}
	public void setStatus(byte status) {
		this.status = status;
	}
	public float getFirstPathRssi() {
		return firstPathRssi;
	}
	public void setFirstPathRssi(float firstPathRssi) {
		this.firstPathRssi = firstPathRssi;
	}
	public float getRssi() {
		return rssi;
	}
	public void setRssi(float rssi) {
		this.rssi = rssi;
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
}
