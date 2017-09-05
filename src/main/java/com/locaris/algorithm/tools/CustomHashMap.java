package com.locaris.algorithm.tools;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import com.locaris.algorithm.beans.TDOADataBean;

public class CustomHashMap<K, V> extends HashMap {
	public CustomHashMap() {
		super();
	}

	public CustomHashMap(int initialCapacity) {
		super(initialCapacity);
	}

	public Object clone() {
		Map target = new HashMap();
		for (Iterator keyIt = this.keySet().iterator(); keyIt.hasNext();) {
			Object key = keyIt.next();
			TDOADataBean tdoaDataBean = (TDOADataBean) this.get(key);
			try {
				target.put(key, tdoaDataBean.clone());
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		return target;
	}
}
