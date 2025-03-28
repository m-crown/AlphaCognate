import { useParams } from "react-router-dom";
import { useEffect, useState } from "react";
import AutoComplete from "../components/SearchBox";
import MolStarViewer from "../components/MolStarViewer";
import TransplantTable from "../components/TransplantTable";
import { Grid, Loader, Text } from "@mantine/core";

type TransplantApiResponse = {
  transplant_data: Array<Transplant>;
  structure_data : Structure;
};

type Transplant = {
  structure_name: string;
  id: number;
  date: string;
  success: boolean;
};

//this is used in multiple places - abstract away
type Structure = {
  name: string;
  url: string;
  transplanted: boolean;
};

function StructurePage() {
  const { id } = useParams(); // Extract the ID from the URL
  const [transplants, setTransplants] = useState<Transplant[]>([]);
  const [structure, setStructure] = useState<Structure | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (id) {
      const fetchTransplantData = async () => {
        const url = new URL(
          '/transplants',
          process.env.NODE_ENV === 'production'
            ? 'http://localhost:8000/' //replace with prod url when it exists
            : 'http://localhost:8000/',
        );
        url.searchParams.append('structure_name', id);
        try {
          const response = await fetch(url.href);
          const json = (await response.json()) as TransplantApiResponse;

          setTransplants(json.transplant_data); // Update state with fetched data
          setStructure(json.structure_data);
          console.log(json.structure);
          console.log(json.transplants);
          setLoading(false); // Set loading to false after data is fetched
        } catch (error) {
          console.error(error);
          setLoading(false); // Set loading to false even if an error occurs
        }
      };
      fetchTransplantData();
    }
  }, [id]);

  return (
    <Grid>
      <Grid.Col span={6}>
        {loading ? <Loader /> : <TransplantTable data={transplants} />}
      </Grid.Col>
      <Grid.Col span={6}>
        <MolStarViewer url={structure?.url || ""} />
      </Grid.Col>
    </Grid>
  );
}

export default StructurePage;
