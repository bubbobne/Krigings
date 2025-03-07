/*

 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.geoframe.blogpost.kriging.primarylocation;


// TODO: Auto-generated Javadoc
/**
 * The Class NumberOfStations.
 */
public class NumberOfStations implements Model{

	/** The in num closer stations. */
	int inNumCloserStations;


	/**
	 * Instantiates a new number of stations.
	 *
	 * @param inNumCloserStations the in number closer stations
	 */
	public NumberOfStations(int inNumCloserStations){

		this.inNumCloserStations=inNumCloserStations;

	}


	/* (non-Javadoc)
	 * @see krigings.Model#numberOfStations()
	 */
	@Override
	public int numberOfStations() {

		return inNumCloserStations;
	}



}
