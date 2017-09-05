package com.locaris.algorithm.pretreatment;

import static org.junit.Assert.assertArrayEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

import com.locaris.algorithm.beans.Position;
import com.locaris.algorithm.beans.TDOADataBean;

public class DeviceFilter {
	private final long Max_TS = 1099511627776L;

	/**
	 * 对基站进行筛选，根据基站时间戳变化量。
	 * 
	 * @param oldData
	 * @param newData
	 * @return
	 */
	public void nlosIdentify(HashMap<Long, TDOADataBean> lastLocData, HashMap<Long, TDOADataBean> thisLocData,
			double vpt) {
		/**
		 * 获取长度
		 */
		if ((lastLocData == null) || (thisLocData == null)) {
			return;
		}
		//将本次定位结果保持一致
		Iterator iterator = thisLocData.entrySet().iterator();
		boolean countFlag=false;
		long timeCompare=0;
		long deviceidFirst=0;
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			if(!countFlag){
				countFlag=true;
				deviceidFirst = (Long) entry.getKey();
				timeCompare= ((TDOADataBean) entry.getValue()).getTimeStamp();
			}else{
				TDOADataBean tdoaDataBean=(TDOADataBean) entry.getValue();
				long temp = tdoaDataBean.getTimeStamp();
				long tempLong = temp - timeCompare;// 获取时间戳变化
				
				if(Math.abs(tempLong)>Max_TS/2){
					if(tempLong>0){
						timeCompare+=Max_TS;
						thisLocData.get(deviceidFirst).setTimeStamp(timeCompare);
					}else{
						temp+=Max_TS;
						tdoaDataBean.setTimeStamp(temp);
					}
				}
			}
		}
		
		iterator = lastLocData.entrySet().iterator();
		countFlag=false;
		timeCompare=0;
		deviceidFirst=0;
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			if(!countFlag){
				countFlag=true;
				deviceidFirst = (Long) entry.getKey();
				timeCompare= ((TDOADataBean) entry.getValue()).getTimeStamp();
			}else{
				TDOADataBean tdoaDataBean=(TDOADataBean) entry.getValue();
				long temp = tdoaDataBean.getTimeStamp();
				long tempLong = temp - timeCompare;// 获取时间戳变化
				if(Math.abs(tempLong)>Max_TS/2){
					if(tempLong>0){
						timeCompare+=Max_TS;
						lastLocData.get(deviceidFirst).setTimeStamp(timeCompare);
					}else{
						temp+=Max_TS;
						tdoaDataBean.setTimeStamp(temp);
					}
				}
			}
		}
		
		
		HashMap<Long, Long> locDataSame = new HashMap<Long, Long>();
		iterator = thisLocData.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			long deviceID = (Long) entry.getKey();
			if (lastLocData.containsKey(deviceID)) {
				TDOADataBean tdoaDataBean = (TDOADataBean) entry.getValue();
				tdoaDataBean.setLastUse(true);
				long thisTimestamp = tdoaDataBean.getTimeStamp();
				long tempLong = thisTimestamp - lastLocData.get(deviceID).getTimeStamp();// 获取时间戳变化
				if (Math.abs(tempLong) > Max_TS / 2) {
					if (tempLong > 0) {
						tempLong = tempLong - Max_TS;
					} else {
						tempLong = tempLong + Max_TS;
					}
				}
				locDataSame.put(deviceID, tempLong);
			}
		}
		if (locDataSame.size() <= 2) {
			return;
		}
		/**
		 * 求取中值
		 */
		long centerdata = 0;
		for (int count = 0; count < locDataSame.size() / 2; count++) {
			long temp = -1;
			iterator = locDataSame.entrySet().iterator();
			while (iterator.hasNext()) {
				Entry entry = (Entry) iterator.next();
				Long timestamp = (Long) entry.getValue();
				if (timestamp > centerdata) {
					if (temp == -1) {
						temp = timestamp;
					} else if (temp > timestamp) {
						temp = timestamp;
					}
				}
			}
			centerdata = temp;
		}
		iterator = locDataSame.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			Long timestamp = (Long) entry.getValue();
			if (Math.abs(timestamp - centerdata) < vpt) {
				// 数据合格
			} else {
				// 数据不合理,将其剔除
				Long deviceID = (Long) entry.getKey();
				thisLocData.remove(deviceID);
			}
		}
	}

	/**
	 * 根据原始数据获取距离本标签最近的基站s。
	 * 
	 * @param thisLocData
	 * @param vpt
	 * @return
	 */
	public HashSet<Long> selectReferencePoints(HashMap<Long, TDOADataBean> thisLocData,
			HashMap<Long, HashSet<Long>> deviceTags, double vpt) {
		HashSet<Long> hashSet = new HashSet<Long>();
		HashSet<Long> nearDevice = new HashSet<Long>();
		long mindata = -1;
		Iterator<Entry<Long, TDOADataBean>> iterator = thisLocData.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry<Long, TDOADataBean> entry = iterator.next();
			TDOADataBean tdoaDataBean = entry.getValue();
			long timestamp = tdoaDataBean.getTimeStamp();
			if (mindata == -1) {
				mindata = timestamp;
			} else if (mindata > timestamp) {
				mindata = timestamp;
			}
		}
		/**
		 * 获取距离合适的基站。
		 */
		iterator = thisLocData.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry<Long, TDOADataBean> entry = iterator.next();
			TDOADataBean tdoaDataBean = entry.getValue();
			if ((tdoaDataBean.getTimeStamp() - mindata) <= vpt) {
				nearDevice.add(entry.getKey());
				// System.out.println("筛选基站ID："+Long.toHexString(entry.getKey()));
				// System.out.println("对应虚拟参考点ID"+deviceTags.get(entry.getKey()));
			}
		}
		/**
		 * 返回所有参考点。
		 */

		for (Iterator it = nearDevice.iterator(); it.hasNext();) {
			HashSet<Long> h = deviceTags.get(it.next());
			if (h != null) {
				hashSet.addAll(h);
			}
		}

		return hashSet;
	}

	/**
	 * 基站筛选----->根据基站间距离以及基站间角度进行基站过滤
	 * 
	 * @param filterData
	 * @param length
	 * @param number
	 * @param distanceVpt
	 * @param angleVpt
	 * @param locInfo
	 * @return
	 */// usestamp 0 id 1 时间戳
	public int lookupFilter(double[][] filterData, long[] deviceids,int length, int number, double distanceVpt, double angleVpt,
			double[][] locInfo,long[] ids) {

		if (length < 3) {
			return 1;// 基站数量不够
		}
		double[][] data = new double[length][5];
		long[] deviceidsc=new long[length];
		for (int i = 0; i < length; i++) {
			System.arraycopy(filterData[i], 0, data[i], 0, 5);
			deviceidsc[i]=deviceids[i];
		}

		// 基站排序,从小到大排序
		for (int i = 0; i < length; i++) {
			for (int j = i + 1; j < length; j++) {
				if (data[i][1] > data[j][1]) {
					double[] temp = data[j];
					data[j] = data[i];
					data[i] = temp;
					long templ=deviceidsc[j];
					deviceidsc[j]=deviceidsc[i];
					deviceidsc[i]=templ;
				}
			}
		}
		/// data为真正使用数据，第一列时间戳，第2.3.4列坐标
		/// 以下代码功能为检索出一组基站距离比较大的基站。
		boolean errflag = false;
		double[][] datastage1 = new double[length][5];
		long[] datastage1ids=new long[length];
		int k = 1;
		// System.arraycopy(data[0], 0, datastage1[0], 0, 4);
		datastage1[0] = data[0];
		datastage1ids[0]=deviceidsc[0];
		for (int i = 1; i < length; i++) {
			for (int j = 0; j < k; j++) {
				double distance = Math.sqrt((data[i][2] - datastage1[j][2]) * (data[i][2] - datastage1[j][2])
						+ (data[i][3] - datastage1[j][3]) * (data[i][3] - datastage1[j][3]));
				if (distance < distanceVpt) {// 如果距离过短
					errflag = true;
					break;
				}
			}
			// 如果所有基站都相距比较远
			if (!errflag) {// 数据可用
				// System.arraycopy(data[i], 0, datastage1[k], 0, 4);
				datastage1[k] = data[i];
				datastage1ids[k]=deviceidsc[i];
				k++;
			}
		}
		if (k < number) {
			// 如果可用基站数量不够
			return 1;// 基站数量不够
		}

		//从大到小排序
		for (int i = 0; i < number; i++) {
			for (int j = i + 1; j < number; j++) {
				if (datastage1[i][1] < datastage1[j][1]) {
					double[] temp = datastage1[j];
					datastage1[j] = datastage1[i];
					datastage1[i] = temp;
					long templ=datastage1ids[j];
					datastage1ids[j]=datastage1ids[i];
					datastage1ids[i]=templ;
				}
			}
		}
		for (int i = 0; i < number; i++) {
			System.arraycopy(datastage1[i], 0, locInfo[i], 0, 5);
			ids[i]=datastage1ids[i];
		}
		return 0;
	}

	/**
	 * 
	 * @param filterData
	 *            ID time x y x
	 * @param deviceNeighbor
	 *            设备ID ，相邻基站
	 * @param locInfo
	 *            ID stamp x y z
	 * @return
	 */
	public int oneDimFilter(long[] deviceids, double[][] filterData, int length,
			HashMap<Long, ArrayList<Long>> deviceNeighbor, double[][] locInfo, long[] ids) {

		double threScale = 0.02;// 比例閾值
		// 时间戳数据小于2，不满足筛选条件.
		if (length < 2) {
			return 0;
		}
		// 找到最小时间戳的基站.
		int minIndex = 0;
		// System.out.println(filterData[length-1][0]);
		int N = length;
		// 获取最小
		for (int i = 0; i < N; i++) {
			if ((filterData[minIndex][0] > filterData[i][0])) {
				minIndex = i;
			}
		}
		long hostId = (deviceids[minIndex]);

		if (!deviceNeighbor.containsKey(hostId)) {
			// System.out.println("hostId");
			return -1;
		}
		double[][] locInfoRaw = new double[3][];// 存储邻域基站定位信息.
		long[] deviceidRaw = new long[3];
		deviceidRaw[0] = deviceids[minIndex];
		locInfoRaw[0] = filterData[minIndex];// 将最近基站对应的一帧数据存储到locInfo.
		int lenSend = 1;
		ArrayList<Long> devices = deviceNeighbor.get(hostId);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < devices.size(); j++) {
				if (deviceids[i] == devices.get(j)) {
					deviceidRaw[lenSend] = deviceids[i];
					locInfoRaw[lenSend] = filterData[i];
					lenSend++;
				}
			}
		}
		if (lenSend < 2) {
			return -2;
		}
		int lenRight = 3;
		boolean[] right = new boolean[lenRight];

		int count = 0;
		for (int i = 0; i < lenSend - 1; i++) {
			for (int j = i + 1; j < lenSend; j++) {
				double tempxx = (locInfoRaw[i][1] - locInfoRaw[j][1]) * (locInfoRaw[i][1] - locInfoRaw[j][1]);
				double tempyy = (locInfoRaw[i][2] - locInfoRaw[j][2]) * (locInfoRaw[i][2] - locInfoRaw[j][2]);
				double tempzz = (locInfoRaw[i][3] - locInfoRaw[j][3]) * (locInfoRaw[i][3] - locInfoRaw[j][3]);

				double dis = Math.sqrt(tempxx + tempyy + tempzz);
				double rdif = (locInfoRaw[i][0] - locInfoRaw[j][0]) / 63897600000.0 * 299792458.0;
				double d1 = (dis + rdif) / 2;
				double d2 = dis - d1;
				double d;
				if (d1 < d2) {
					d = d1;
				} else {
					d = d2;
				}
				double scale = d / dis;
				if (scale > threScale) {
					right[count] = true;
				} else {
					right[count] = false;
				}
				count = count + 1;
			}
		}
		if (count == 1) {
			if (right[0] == true) {
				locInfo[0] = locInfoRaw[0];
				locInfo[1] = locInfoRaw[1];

				return 2;
			} else {
//				 System.out.println("count11111");
				return -3;
			}

		}
		int id1;
		int id2;
		int id3;
		if (count == 3) {
			if (right[0] == true && right[1] == false) {
				id1 = 0;
				id2 = 1;
				id3 = 2;
			} else if (right[0] == false && right[1] == true) {
				id1 = 0;
				id2 = 2;
				id3 = 1;
			} else if (right[0] == false && right[1] == false && right[2] == true) {
				id1 = 1;
				id2 = 2;
				id3 = 0;
			} else if (right[0] == true && right[1] == true) {
				id1 = 1;
				id2 = 2;
				id3 = 0;
			} else {
				// System.out.println("count22222");
				// System.out.println(right[0]);
				// System.out.println(right[1]);
				// System.out.println(right[2]);
				return -4;
			}
			locInfo[0] = locInfoRaw[id1];
			locInfo[1] = locInfoRaw[id2];
			locInfo[2] = locInfoRaw[id3];
			ids[0] = deviceidRaw[id1];
			ids[1] = deviceidRaw[id2];
			ids[2] = deviceidRaw[id3];
			return 3;
		}
		// System.out.println("count33333");
		return -5;
	}

}
