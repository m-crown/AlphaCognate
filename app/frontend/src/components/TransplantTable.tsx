import { useMemo, useState } from 'react';
import {
  MantineReactTable,
  useMantineReactTable,
  type MRT_ColumnDef,
  type MRT_ColumnFiltersState,
  type MRT_PaginationState,
  type MRT_SortingState,
} from 'mantine-react-table';

type Transplant = {
  id: string;
  structure_name: string;
  date: string;
};

interface TransplantTableProps {
  data: Transplant[]; // Accept data as a prop
}

const TransplantTable = ({ data }: TransplantTableProps) => {
  console.log("transplants")
  console.log(data) 
  // table state
  const [columnFilters, setColumnFilters] = useState<MRT_ColumnFiltersState>([]);
  const [globalFilter, setGlobalFilter] = useState('');
  const [sorting, setSorting] = useState<MRT_SortingState>([]);
  const [pagination, setPagination] = useState<MRT_PaginationState>({
    pageIndex: 0,
    pageSize: 10,
  });

  const columns = useMemo<MRT_ColumnDef<Transplant>[]>(
    () => [
      {
        accessorKey: 'id',
        header: 'Transplant ID',
      },
      {
        accessorKey: 'structure_name',
        header: 'Structure Name',
      },
      {
        accessorKey: 'date',
        header: 'Date',
      },
    ],
    [],
  );

  const table = useMantineReactTable({
    columns,
    data,
    enableRowSelection: true,
    getRowId: (row) => row.id,
    initialState: { showColumnFilters: true },
    manualFiltering: true,
    manualPagination: true,
    manualSorting: true,
    onColumnFiltersChange: setColumnFilters,
    onGlobalFilterChange: setGlobalFilter,
    onPaginationChange: setPagination,
    onSortingChange: setSorting,
    state: {
      columnFilters,
      globalFilter,
      pagination,
      sorting,
    },
  });

  return (
    <div style={{ maxHeight: '80vh', maxWidth: '100%', overflow: 'auto' }}>
      <MantineReactTable table={table} />
    </div>
  );
};

export default TransplantTable;
