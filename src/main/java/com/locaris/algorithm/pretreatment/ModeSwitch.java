package com.locaris.algorithm.pretreatment;

import static org.hamcrest.CoreMatchers.instanceOf;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

import com.locaris.algorithm.beans.DevicePos;
import com.locaris.algorithm.beans.TDOADataBean;
import com.locaris.algorithm.beans.TagPos;
import com.locaris.algorithm.beans.TagResidual;

/**
 * 定位维度切换模块
 * 
 * @author Lyn
 *
 */
public class ModeSwitch {
	/**
	 * 模式切换初始化,获取特征值
	 * @param deviceHashMap--->基站坐标嘻哈表
	 * @param tagHashMap------>标签坐标嘻哈表
	 * @param tagDistance----->结果,距离
	 */
	public void switchInit(HashMap<Long, DevicePos> deviceHashMap, HashMap<Long, TagPos> tagHashMap,
			HashMap<Long, HashMap<Long, Double>> tagDistance,HashMap<Long, HashSet<Long>> deviceTags) {
		Iterator iterator = tagHashMap.entrySet().iterator();
		tagDistance.clear();
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			Long tagID = (Long) entry.getKey();
			TagPos tagPos = (TagPos) entry.getValue();
			Iterator iterator2 = deviceHashMap.entrySet().iterator();
			double minLength=-1;
			long minid=-1;
			HashMap<Long, Double> hashMap=new HashMap<Long, Double>();
			while (iterator2.hasNext()) {
				Entry entry2 = (Entry) iterator2.next();
				DevicePos devicePos = (DevicePos) entry2.getValue();
				Long deviceID = (Long) entry2.getKey();
				double temp = (tagPos.getX() - devicePos.getX()) * (tagPos.getX() - devicePos.getX())
						+ (tagPos.getY() - devicePos.getY()) * (tagPos.getY() - devicePos.getY())
						+ (tagPos.getZ() - devicePos.getZ()) * (tagPos.getZ() - devicePos.getZ());
				if(minLength==-1){
					minLength=temp;
					minid=deviceID;
				}else if(minLength>temp){
					minLength=temp;
					minid=deviceID;
				}
				hashMap.put(deviceID, Math.sqrt(temp));	
			}
			/**
			 * 新添加，找到device对应的最小距离的参考点。
			 */
//			System.out.println("筛选基站ID："+Long.toHexString(minid));
//			System.out.println("对应虚拟参考点ID："+Long.toHexString(tagID));
			if(minid!=-1){
				if(deviceTags.containsKey(minid)){
					HashSet<Long> hashSet=deviceTags.get(minid);
					hashSet.add(tagID);
				}else{
					HashSet<Long> hashSet=new HashSet<Long>();
					hashSet.add(tagID);
					deviceTags.put(minid, hashSet);
				}
			}
			tagDistance.put(tagID, hashMap);
			Iterator iterator3 = hashMap.entrySet().iterator();
			while (iterator3.hasNext()) {
				Entry entry3 = (Entry) iterator3.next();
				entry3.setValue((Double)entry3.getValue()-Math.sqrt(minLength));
			}
		}
	}

	/**
	 * 切换定位模式
	 * @param locInfo-->定位信息 第一列:基站ID,第二列:对应时间戳     时间戳单位统一为米
	 * @param tagDistance
	 * @param length
	 * @param minLength-->最少需要几个基站的数据
	 * @return
	 */
	public TagResidual switchMode(HashMap<Long,TDOADataBean> locInfo,HashMap<Long, HashMap<Long, Double>> tagDistance,HashSet<Long> hashSet,int minLength) {
		int length=locInfo.size();
		if(length<minLength){
			return null;
		}
		if(length>3){
			length=3;
		}
		double[][] locInfoCopy=new double[length][2];
		/**
		 * 获取最小时间戳,并求取距差
		 */
		long minTimestamp=-1;
		Iterator iterator = locInfo.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			Long timestamp = ((TDOADataBean) entry.getValue()).getTimeStamp();
			if(minTimestamp==-1){
				minTimestamp=timestamp;
			}else if (minTimestamp > timestamp) {
				minTimestamp=timestamp;
			}
		}
		iterator = locInfo.entrySet().iterator();
		int icount=0;
		while (iterator.hasNext()) {
			Entry entry = (Entry) iterator.next();
			Long timestamp = ((TDOADataBean) entry.getValue()).getTimeStamp();
			locInfoCopy[icount][0]=(Long) entry.getKey();
			locInfoCopy[icount][1]=timestamp-minTimestamp;
			locInfoCopy[icount][1] = locInfoCopy[icount][1] / 63897600000L * 299792458L;
			icount++;
			if(icount==length){
				break;
			}
		}
		/**
		 * 根据距差从小到大排序
		 */
		for(int i=0;i<length;i++){
			for(int j=i+1;j<length;j++){
				if (locInfoCopy[i][1] > locInfoCopy[j][1]) {
					double[] temp=locInfoCopy[i];
					locInfoCopy[i]=locInfoCopy[j];
					locInfoCopy[j]=temp;
	            }
			}
		}
		double minResidual=-1;
		Long minTagID=(long) -1;
		/**
		 * 遍历参考标签点,并求其残差.
		 */
		for (Iterator it = hashSet.iterator(); it.hasNext();) {
			double residual=0;
			long referenctID=(Long) it.next();
			HashMap<Long, Double> hashMap=tagDistance.get(referenctID);
			for(int count=0;count<length;count++){
				Long deviceID=(long) locInfoCopy[count][0];
				residual+=(locInfoCopy[count][1]-hashMap.get(deviceID))*(locInfoCopy[count][1]-hashMap.get(deviceID));
			}
			residual/=length;
			/**
			 * 最小残差
			 */
			if(minResidual==-1){
				minResidual=residual;
				minTagID=referenctID;
			}else if(minResidual>residual){
				minResidual=residual;
				minTagID=referenctID;
			}
		}
		if((minResidual==-1)||(minTagID==-1)){
			return null;
		}
		//返回最小的参考标签残差值与其对应的标签ID
		TagResidual tagResidual=new TagResidual();
		tagResidual.setTagID(minTagID);
		tagResidual.setResidual(minResidual);
		return tagResidual;
	}
}
