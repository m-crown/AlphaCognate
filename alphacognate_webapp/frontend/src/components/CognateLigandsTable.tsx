import { useMemo } from 'react';
import {
  MantineReactTable,
  type MRT_ColumnDef
} from 'mantine-react-table';

interface CognateLigand {
  id: number;
  name: string;
  smiles: string;
  xref: string;
  similarity: number;
}

interface CognateLigandsTableProps {
  cognateLigandsData: CognateLigand[];
}

const CognateLigandsTable = ({ cognateLigandsData }: CognateLigandsTableProps) => {
  const columns = useMemo<MRT_ColumnDef<CognateLigand>[]>(
    () => [
      {
        accessorKey: 'name',
        header: 'Cognate Ligand Name',
      },
      {
        accessorKey: 'smiles',
        header: 'SMILES',
      },
      {
        accessorKey: 'xref',
        header: 'Cross References',
      },
      {
        accessorKey: 'similarity',
        header: 'Similarity',
      },
    ],
    [],
  );

  return (
    <div style={{ maxHeight: '80vh', overflow: 'auto' }}>
      <MantineReactTable
        columns={columns}
        data={cognateLigandsData}
        initialState={{ showColumnFilters: true }}
      />
    </div>
  );
};

export default CognateLigandsTable;
