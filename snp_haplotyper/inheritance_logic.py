from abc import ABC, abstractmethod
from io import StringIO
from typing import Dict, List, Tuple, Union

import pandas as pd
from EnumDataClasses import InheritanceMode, Relationship, Sex, Status
from helper_functions import set_inherited_from_category_dtype, set_risk_category_dtype


class InheritanceLogic(ABC):
    def __init__(
        self,
        df: pd.DataFrame,
        partner1: str,
        partner2: str,
        reference: str,
        reference_status: Status,
        reference_relationship: Relationship,
        reference_sex: Sex,
        consanguineous: bool,
    ):
        self.df = df
        self.partner1 = partner1  # AD = affected_partner, AR = female_partner, XL = carrier_female_partner
        self.partner2 = partner2  # AD = unaffected_partner, AR = male_partner, XL = unaffected_male_partner
        self.reference = reference
        self.reference_status = reference_status
        self.reference_relationship = reference_relationship
        self.reference_sex = reference_sex
        self.consanguineous = consanguineous
        self.lookup_df: pd.DataFrame = self.select_lookup_df()
        self.add_nocall_lookup()
        self.merge_dataframes()

    def add_nocall_lookup(self):
        self.nocall_data = f"""
        {self.reference}, {self.partner1},  {self.partner2},    snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               AB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               BB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               NoCall,           AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               NoCall,           AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               NoCall,           BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AA,               NoCall,           NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               AA,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               AB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               BB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               NoCall,           AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               NoCall,           AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               NoCall,           BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        AB,               NoCall,           NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               AA,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               AB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               BB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               NoCall,           AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               NoCall,           AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               NoCall,           BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        BB,               NoCall,           NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AA,               AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AA,               AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AA,               BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AA,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AB,               AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AB,               AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AB,               BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           AB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           BB,               AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           BB,               AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           BB,               BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           BB,               NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           NoCall,           AA,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           NoCall,           AB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           NoCall,           BB,                 NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        NoCall,           NoCall,           NoCall,             NoCall_in_trio,       NoCall_in_trio,       NoCall_in_trio,       unassigned
        """
        if self.nocall_data is None:
            raise ValueError("Data not provided in the subclass")

        data_io = StringIO(self.nocall_data)
        try:
            # Assuming the CSV data is well-formed and includes headers
            nocall_lookup_df = pd.read_csv(
                data_io, sep=",", comment="#", skipinitialspace=True
            )
            nocall_lookup_df["snp_inherited_from"] = nocall_lookup_df[
                "snp_inherited_from"
            ].str.strip()  # Remove whitespace from snp_inherited_from column due to Comments
            self.lookup_df = pd.concat([self.lookup_df, nocall_lookup_df])
        except Exception as e:
            # Handle any exceptions that might occur while reading the CSV data
            raise ValueError(f"Error reading NoCalldata: {e}")
        # Append to the lookup_df

    def merge_dataframes(self):
        # Merging self.df with self.lookup_df
        merge_cols = [self.reference, self.partner1, self.partner2]
        self.df = pd.merge(
            self.df, self.lookup_df, on=merge_cols, how="left", suffixes=("", "_lookup")
        )

        self.df = set_risk_category_dtype(self.df)
        self.df = set_inherited_from_category_dtype(self.df)

    @abstractmethod
    def select_lookup_df(self):
        lookup_df: pd.DataFrame = pd.DataFrame()
        return lookup_df

    def get_results(self):
        return self.df


class AutosomalDominantLogic(InheritanceLogic):
    def __init__(
        self,
        df: pd.DataFrame,
        affected_partner: str,
        unaffected_partner: str,
        reference: str,
        reference_status: Status,
        reference_relationship: Relationship,
        reference_sex: Sex,
        consanguineous: bool,
    ):
        super().__init__(
            df,
            affected_partner,
            unaffected_partner,
            reference,
            reference_status,
            reference_relationship,
            reference_sex,
            consanguineous,
        )
        # Additional initialization for AutosomalDominantLogic
        self.df = df
        self.partner1 = affected_partner  # AD = affected_partner, AR = female_partner, XL = carrier_female_partner
        self.partner2 = unaffected_partner  # AD = unaffected_partner, AR = male_partner, XL = unaffected_male_partner
        self.reference = reference
        self.reference_status = reference_status
        self.reference_relationship = reference_relationship
        self.reference_sex = reference_sex
        self.consanguineous = consanguineous
        self.look_up_df: pd.DataFrame = self.select_lookup_df()
        self.add_nocall_lookup()
        self.merge_dataframes()

    def select_lookup_df(self):
        lookup_df = pd.DataFrame()
        if self.reference_status == Status.AFFECTED:
            if self.reference_relationship == Relationship.GRANDPARENT:
                lookup_df = AD_RefAffectedGrandparentLookup().get_dataframe()
            elif self.reference_relationship == Relationship.CHILD:
                lookup_df = AD_RefAffectedChildLookup().get_dataframe()
        elif self.reference_status == Status.UNAFFECTED:
            if self.reference_relationship == Relationship.GRANDPARENT:
                lookup_df = AD_RefUnaffectedGrandparentLookup().get_dataframe()
            elif self.reference_relationship == Relationship.CHILD:
                lookup_df = AD_RefUnaffectedChildLookup().get_dataframe()
        # if lookup_df is not defined, raise an error
        if lookup_df.empty:
            raise ValueError(
                "lookup_df is empty, check reference_status and reference_relationship"
            )
        # Rename the columns to match the input dataframe
        lookup_df.rename(
            columns={
                "reference": self.reference,
                "affected_partner": self.partner1,
                "unaffected_partner": self.partner2,
            },
            inplace=True,
        )

        return lookup_df


class AutosomalRecessiveLogic(InheritanceLogic):
    def __init__(
        self,
        df: pd.DataFrame,
        male_partner: str,
        female_partner: str,
        reference: str,
        reference_status: Status,
        reference_relationship: Relationship,
        reference_sex: Sex,
        consanguineous: bool,
    ):
        super().__init__(
            df,
            male_partner,
            female_partner,
            reference,
            reference_status,
            reference_relationship,
            reference_sex,
            consanguineous,
        )

        # Additional initialization for AutosomalRecessiveLogic
        self.df = df
        self.partner1 = female_partner  # AD = affected_partner, AR = female_partner, XL = carrier_female_partner
        self.partner2 = male_partner  # AD = unaffected_partner, AR = male_partner, XL = unaffected_male_partner
        self.reference = reference
        self.reference_status = reference_status
        self.reference_relationship = reference_relationship
        self.reference_sex = reference_sex
        self.consanguineous = consanguineous
        self.look_up_df: pd.DataFrame = self.select_lookup_df()
        self.add_nocall_lookup()
        self.merge_dataframes()

    def select_lookup_df(self):
        lookup_df = pd.DataFrame()
        if self.reference_status == Status.AFFECTED:
            if self.consanguineous:
                lookup_df = AR_RefAffectedConsangLookup().get_dataframe()
            else:
                lookup_df = AR_RefAffectedNonconsangLookup().get_dataframe()
        elif self.reference_status == Status.UNAFFECTED:
            if self.consanguineous:
                print("Cannot have unaffected reference in consanguineous analysis")
            else:
                lookup_df = AR_RefUnaffectedNonconsangLookup().get_dataframe()
        # Rename the columns to match the input dataframe
        lookup_df.rename(
            columns={
                "reference": self.reference,
                "male_partner": self.partner2,
                "female_partner": self.partner1,
            },
            inplace=True,
        )
        return lookup_df


class XLinkedLogic(InheritanceLogic):
    def __init__(
        self,
        df: pd.DataFrame,
        carrier_female_partner: str,
        unaffected_male_partner: str,
        reference: str,
        reference_status: Status,
        reference_relationship: Relationship,
        reference_sex: Sex,
        consanguineous: bool,
    ):
        super().__init__(
            df,
            carrier_female_partner,
            unaffected_male_partner,
            reference,
            reference_status,
            reference_relationship,
            reference_sex,
            consanguineous,
        )
        # Additional initialization for XLinkedLogic

        self.df = df
        self.partner1 = carrier_female_partner  # AD = affected_partner, AR = female_partner, XL = carrier_female_partner
        self.partner2 = unaffected_male_partner  # AD = unaffected_partner, AR = male_partner, XL = unaffected_male_partner
        self.reference = reference
        self.reference_status = reference_status
        self.reference_relationship = reference_relationship
        self.reference_sex = reference_sex
        self.consanguineous = consanguineous
        self.look_up_df: pd.DataFrame = self.select_lookup_df()
        self.add_nocall_lookup()
        self.merge_dataframes()

    def select_lookup_df(self):
        lookup_df = pd.DataFrame()
        if self.reference_sex == Sex.MALE:
            lookup_df = XL_RefMaleLookup().get_dataframe()
        elif self.reference_sex == Sex.FEMALE:
            lookup_df = XL_RefFemaleLookup().get_dataframe()
        elif self.reference_sex == Sex.UNKNOWN:
            print("X-Linked References must be sexed")
        else:
            print("Sex must be Sex.MALE, Sex.FEMALE, or Sex.UNKNOWN")
        # Rename the columns to match the input dataframe
        lookup_df.rename(
            columns={
                "reference": self.reference,
                "carrier_female_partner": self.partner1,
                "unaffected_male_partner": self.partner2,
            },
            inplace=True,
        )

        return lookup_df


class BaseLookup(ABC):
    @abstractmethod
    def __init__(self):
        # Subclasses should initialize their own data
        self.data = None

    def get_dataframe(self):
        """
        Converts the data string to a pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame created from the data string.
        """
        # Check if data is provided by the subclass
        if self.data is None:
            raise ValueError("Data not provided in the subclass")

        data_io = StringIO(self.data)
        try:
            # Assuming the CSV data is well-formed and includes headers
            df = pd.read_csv(data_io, sep=",", comment="#", skipinitialspace=True)
            df["snp_inherited_from"] = df[
                "snp_inherited_from"
            ].str.strip()  # Remove whitespace from snp_inherited_from column due to Comments
            return df
        except Exception as e:
            # Handle any exceptions that might occur while reading the CSV data
            raise ValueError(f"Error reading data: {e}")


class AD_RefAffectedGrandparentLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Dominant
        # Reference_status: affected
        # Reference_relationship: grandparent
        self.data = """
        reference,        affected_partner, unaffected_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 7
        AA,               AB,               BB,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,               AA,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 25
        AB,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AD_RefUnaffectedGrandparentLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Dominant
        # Reference_status: unaffected
        # Reference_relationship: grandparent
        self.data = """
        reference,        affected_partner, unaffected_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 7
        AA,               AB,               BB,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,               AA,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 25
        AB,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AD_RefAffectedChildLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Dominant
        # Reference_status: affected
        # Reference_relationship: child or embryo
        self.data = """
        reference,        affected_partner, unaffected_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 7
        AA,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,               AA,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 25
        AB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AD_RefUnaffectedChildLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Dominant
        # Reference_status: unaffected
        # Reference_relationship: child or embryo
        self.data = """
        reference,        affected_partner, unaffected_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 7
        AA,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        unassigned         # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,               AA,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 25
        AB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        unassigned         # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AR_RefAffectedNonconsangLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Recessive
        # Reference_status: affected
        # Consanguineous: False
        self.data = """
        reference,        male_partner,     female_partner,       snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 7
        AA,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 15
        BB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 24
        AB,               AB,               AA,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 25
        AB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AR_RefUnaffectedNonconsangLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Recessive
        # Reference_status: affected
        # Consanguineous: False
        self.data = """
        reference,        male_partner,     female_partner,     snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 7
        AA,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 15
        BB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 17
        BB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 24
        AB,               AB,               AA,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 25
        AB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class AR_RefAffectedConsangLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: Autosomal Recessive
        # Reference_status: affected
        # Consanguineous: True
        self.data = """
        reference,        male_partner,     female_partner,     snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 3
        AA,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,               AA,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 7
        AA,               AB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 8
        AA,               AB,               AB,                 uninformative,        high_risk,            low_risk,             both_partners      # Comment for row 9
        BB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,               AB,                 low_risk,             uninformative,        uninformative,        female_partner     # Comment for row 15
        BB,               AB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 16
        BB,               AB,               BB,                 low_risk,             uninformative,        uninformative,        male_partner       # Comment for row 17
        BB,               AB,               AB,                 uninformative,        low_risk,             high_risk,            both_partners      # Comment for row 18
        AB,               AA,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 21
        AB,               BB,               AA,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,               BB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,               AB,                 high_risk,            uninformative,        uninformative,        female_partner     # Comment for row 24
        AB,               AB,               AA,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 25
        AB,               AB,               BB,                 high_risk,            uninformative,        uninformative,        male_partner       # Comment for row 26
        AB,               AB,               AB,                 uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class XL_RefMaleLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: X-linked
        # Reference_sex: Male
        self.data = """
        reference,        carrier_female_partner, unaffected_male_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,                     AA,                   low_risk,             high_risk,            low_risk,             unassigned         # Comment for row 7
        AA,               AB,                     BB,                   high_risk,            high_risk,            low_risk,             unassigned         # Comment for row 8
        AA,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,                     AA,                   high_risk,            low_risk,             high_risk,            unassigned         # Comment for row 16
        BB,               AB,                     BB,                   low_risk,             low_risk,             high_risk,            unassigned         # Comment for row 17
        BB,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 25
        AB,               AB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 26
        AB,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """


class XL_RefFemaleLookup(BaseLookup):
    def __init__(self):
        super().__init__()
        # InheritanceMode: X-linked
        # Reference_sex: Female
        self.data = """
        reference,        carrier_female_partner, unaffected_male_partner, snp_risk_category_AB, snp_risk_category_AA, snp_risk_category_BB, snp_inherited_from
        AA,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 1
        AA,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 2
        AA,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 3
        AA,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 4
        AA,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 5
        AA,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 6
        AA,               AB,                     AA,                   low_risk,             high_risk,            low_risk,             unassigned         # Comment for row 7
        AA,               AB,                     BB,                   high_risk,            high_risk,            low_risk,             unassigned         # Comment for row 8
        AA,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 9
        BB,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 10
        BB,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 11
        BB,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 12
        BB,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 13
        BB,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 14
        BB,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 15
        BB,               AB,                     AA,                   high_risk,            low_risk,             high_risk,            unassigned         # Comment for row 16
        BB,               AB,                     BB,                   low_risk,             low_risk,             high_risk,            unassigned         # Comment for row 17
        BB,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 18
        AB,               AA,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 19
        AB,               AA,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 20
        AB,               AA,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 21
        AB,               BB,                     AA,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 22
        AB,               BB,                     BB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 23
        AB,               BB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 24
        AB,               AB,                     AA,                   high_risk,            low_risk,             high_risk,            unassigned         # Comment for row 25
        AB,               AB,                     BB,                   high_risk,            high_risk,            low_risk,             unassigned         # Comment for row 26
        AB,               AB,                     AB,                   uninformative,        uninformative,        uninformative,        unassigned         # Comment for row 27
        """
